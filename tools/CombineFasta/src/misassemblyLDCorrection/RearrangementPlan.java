/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package misassemblyLDCorrection;

import de.erichseifert.gral.data.DataSeries;
import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.graphics.Drawable;
import de.erichseifert.gral.graphics.DrawingContext;
import de.erichseifert.gral.graphics.Label;
import de.erichseifert.gral.plots.XYPlot;
import de.erichseifert.gral.plots.colors.LinearGradient;
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.lines.LineRenderer;
import de.erichseifert.gral.plots.points.DefaultPointRenderer2D;
import de.erichseifert.gral.plots.points.PointRenderer;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import javax.imageio.ImageIO;

/**
 * The main workhorse for identifying which chromosomes and segments need to be rearranged
 * Reads in an alignment file and generates segment "blocks" that are to be reordered
 * @author dbickhart
 */
public class RearrangementPlan {
    private static final Logger log = Logger.getLogger(RearrangementPlan.class.getName());
    private final Path samFile;
    private IndexedFastaReader origin; 
    //private BufferedFastaReaderWriter output;
    //protected static final Pattern readMarkerOrder = Pattern.compile(".+\\.(.+)\\.(.+)");
    private List<markerCoords> coords;
    private Map<String, List<BedFastaPlan>> mergedPlan = new HashMap<>();
    private Set<String> unmodified;
    //private static final MergerUtils util = new MergerUtils();
    
    public RearrangementPlan(String samFile, String fasta){
        this.samFile = Paths.get(samFile);
        if(!this.samFile.toFile().canRead()){
            log.log(Level.SEVERE, "Cannot read input sam file! terminating...");
            System.exit(-1);
        }
        
        this.origin = new IndexedFastaReader(Paths.get(fasta));
    }
    
    public void CreateMarkerPlan(){
        // Read and process the sam file into a list        
        try(BufferedReader input = Files.newBufferedReader(this.samFile, Charset.defaultCharset())){
            this.coords = input.lines()
                    .filter((s) -> (! s.startsWith("@")))
                    .map((l) -> {return new markerCoords(l);})
                    .filter((m) -> !m.isUnmapped())
                    .collect(Collectors.toList());
        }catch(IOException ex){
            log.log(Level.SEVERE, "Encountered error reading sam file!", ex);
        }
        
        log.log(Level.INFO, "Finished sam file reading and organization");
                
        // sort the list and process into a bed file
        // First, sort by new marker maps and link prior markers together
        coords.sort(new NewChrSort()
                .thenComparingInt(markerCoords::getNPos));
        log.log(Level.INFO, "Finished sort by origin fasta coordinates");
        
        // Fill in the linked list for each marker based on new position coordinates
        coords.get(0).setNextMarker(coords.get(1));
        coords.get(0).setOrder(0);
        for(int x = 1; x < coords.size() - 1; x++){
            coords.get(x).setPrevMarker(
                    (coords.get(x - 1).nChr.equals(coords.get(x).nChr))? coords.get(x - 1) : null);
            coords.get(x).setNextMarker(
                    (coords.get(x + 1).nChr.equals(coords.get(x).nChr)? coords.get(x + 1) : null));
            coords.get(x).setOrder(x);
        }
        coords.get(coords.size() - 1).setOrder(coords.size() -1);
        
        // Now sort by original coordinate order
        // Such a better syntax for comparisons here! 
        coords.sort(new OrgChrSort()
                .thenComparingInt(markerCoords::getOPos));
        log.log(Level.INFO, "Finished sort by marker order coordinates");
        
        // Finally, merge into overlapping segments
        String lastChr = coords.get(0).oChr;
        this.mergedPlan.put(lastChr, new ArrayList<>());
        this.mergedPlan.get(lastChr).add(convertCoordToBed(coords.get(0)));
        
        for(int x = 1; x < coords.size(); x++){
            markerCoords current = coords.get(x);
            if(this.mergedPlan.containsKey(current.oChr)){
                BedFastaPlan ref = this.mergedPlan.get(current.oChr)
                        .get(this.mergedPlan.get(current.oChr).size() - 1);
                BedFastaPlan query = convertCoordToBed(current);
                //if(MergerUtils.checkOverlap(ref.Start(), query.Start(), ref.End(), query.Start())
                //        && ref.Chr().equals(query.Chr())){
                if(Math.abs(current.order - coords.get(x - 1).order) == 1
                        && ref.Chr().equals(query.Chr())){
                    // Check if this is a rev comp join
                    if(current.order - coords.get(x - 1).order < 0){
                        // The merger is a rev comp join
                        if(query.Start() > ref.End())
                            log.log(Level.WARNING, "REV logic! Attempted to refine end of: " + ref.toString() + " to " + query.toString());
                        ref.setStart(query.Start());
                        ref.isRev = true;
                    }else{
                        if(query.End() < ref.Start())
                            log.log(Level.WARNING, "BAD logic! Attempted to refine end of: " + ref.toString() + " to " + query.toString());
                        ref.setEnd(query.End());
                    }
                    ref.incCounter();
                }else{
                    // No overlap, so we'll add this as a separate entry
                    this.mergedPlan.get(current.oChr).add(query);
                }
            }else{
                // Need to start a new chromosome
                this.mergedPlan.put(current.oChr, new ArrayList<>());
                this.mergedPlan.get(current.oChr).add(convertCoordToBed(current));
            }
            lastChr = current.oChr;
        }
        log.log(Level.INFO, "Finished marker plan");
        
        // Merge smaller events
        // Temporarily disabled if the conditional is uncommented
        //if(false){
            this.mergedPlan.entrySet().stream().forEach((s) -> {
                String chr = s.getKey();
                List<BedFastaPlan> existing = s.getValue();
                List<BedFastaPlan> merged = new ArrayList<>();
                for(int x = 1; x < existing.size() - 1; x++){
                    BedFastaPlan current = existing.get(x);
                    BedFastaPlan previous = existing.get(x - 1);
                    BedFastaPlan future = existing.get(x + 1);
                    if(current.counter < 2 && previous.counter > 2 && future.counter > 2){
                        // Small segment between two larger ones
                        if(previous.End() + 50000 < future.Start() && previous.Chr().equals(future.Chr())
                                && previous.isRev == future.isRev){
                            int end = utils.MergerUtils.most(previous.End(), future.End());
                            int start = utils.MergerUtils.least(previous.Start(), future.End());
                            BedFastaPlan temp = new BedFastaPlan(previous.Chr(), 
                                    start, end, previous.getActualChr(), previous.getActualStart());
                            temp.isRev = previous.isRev;
                            merged.add(temp);
                            x += 2;
                            continue;
                        }
                    }else if(!previous.Chr().equals(chr) && previous.counter < 5){
                        // Small segment that may be a mismapped marker
                        continue;
                    }
                    merged.add(previous);
                }
                int sizeDiff = existing.size() - merged.size();
                log.log(Level.INFO, "[MERGESEGS] Original chr: " + chr + "\tmerged " + sizeDiff + " smallers segments.");
                this.mergedPlan.put(chr, merged);
            });
        //}
        
        // Generate current chromosome segment stats
        this.mergedPlan.entrySet().stream().forEach((s) -> {
            String chr = s.getKey();
            int num = s.getValue().size();
            long diffChrs = s.getValue().stream()
                    .map(BedFastaPlan::Chr)
                    .count();
            log.log(Level.INFO, "[MERGESTATS] Original chr: " + chr + "\tsegments: " + num + "\tsegment chr counts: " + diffChrs);
        });
        
        // Identify unmodified chromosomes
        this.unmodified = this.origin.getChrNames().stream()
                .filter(s -> !this.mergedPlan.containsKey(s))
                .collect(Collectors.toSet());
        
        log.log(Level.INFO, "[MERGESTATS] Number of chrs unchanged: " + this.unmodified.size());
        
        /*//Ensuring that the list is cleared from memory
        Removed so that a separate plot method could be derived
        this.coords = null;
        System.gc();*/
    }
    
    public void plotGraph(String outfile){
        // Assign data to table and count number of elements
        int numElements = this.coords.size();
        Map<String, Long> chrElementCount = this.coords.stream()
                .map(s -> s.oChr)
                .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()));
        
        // Sorting chromosome segment order by current coord chr set
        List<String> chrOrder = this.coords.stream()
                .map(markerCoords::getOChr)
                .collect(Collectors.toSet())
                .stream().sorted((s, s1) ->{
                            if(utils.MergerUtils.isNumeric(s) 
                                    && utils.MergerUtils.isNumeric(s1)){
                                return Integer.compare(Integer.parseInt(s), Integer.parseInt(s1));
                            }
                            return s.compareTo(s1);}).collect(Collectors.toList());
        
        // Generating table and data series
        DataTable data = new DataTable(Integer.class, Integer.class);
        for(int x = 0; x < numElements; x++){
            data.add(x, this.coords.get(x).order);
        }
        
        DataSeries series1 = new DataSeries("Marker Order", data, 0, 1);
        XYPlot plot = new XYPlot(series1);
        plot.getTitle().setText("Relative marker order in assembly");
        
        
        // Style the points to be a small rectangle
        PointRenderer points = new DefaultPointRenderer2D();
        points.setShape(new Rectangle2D.Double(-2, -2, 2, 2));
        points.setColor(new LinearGradient(Color.BLUE, Color.CYAN, Color.DARK_GRAY, Color.GRAY, Color.lightGray, Color.GREEN, Color.YELLOW, Color.ORANGE, Color.RED));
        plot.setPointRenderers(data, points);
        
        // Chromosome boundary lines
        Map<Double, String> tickPositions = new HashMap<>();
        double prevPos = 0;
        for(String chr : chrOrder){
            tickPositions.put(chrElementCount.get(chr) + prevPos, chr);
            prevPos += chrElementCount.get(chr);
            log.log(Level.INFO, "[PLOTTING] Chr " + chr + " tick mark set at: " + prevPos);
        } 
        
        plot.getAxisRenderer(XYPlot.AXIS_X).setLabel(new Label("Original Marker Order"));
        plot.getAxisRenderer(XYPlot.AXIS_Y).setLabel(new Label("New Marker Order"));
        plot.getAxisRenderer(XYPlot.AXIS_X).setTicksAutoSpaced(false);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setTicksAutoSpaced(false);
        
        plot.getAxisRenderer(XYPlot.AXIS_X).setCustomTicks(tickPositions);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setCustomTicks(tickPositions);
        plot.getAxisRenderer(XYPlot.AXIS_X).setTickLabelsVisible(true);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setTickLabelsVisible(true);
        plot.getAxisRenderer(XYPlot.AXIS_X).setTickLabelsOutside(false);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setTickLabelsOutside(false); 
        plot.getAxisRenderer(XYPlot.AXIS_X).setIntersection(-Double.MAX_VALUE);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setIntersection(-Double.MAX_VALUE);
        plot.getAxisRenderer(XYPlot.AXIS_X).setMinorTicksVisible(false);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setMinorTicksVisible(false);
        
        writePNG(plot, numElements, new File(outfile));
        
        // After writing, we can clear the original marker elements
        this.coords = null;
        System.gc();
    }
    
    private void writePNG(Drawable data, int numEntries, File outputFile){
        int width = 8000, height = 8000;
        
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D graphics = image.createGraphics();
        graphics.setPaintMode();
        
        DrawingContext ctx = new DrawingContext(graphics);
        data.setBounds(0, 0, width, height);
        data.draw(ctx);
        
        try {
            ImageIO.write(image, "png", outputFile);
        } catch (IOException ex) {
            log.log(Level.SEVERE, "Error writing plot to PNG!", ex);
        }
    }
    
    public void refineEdges(String JellyfishDb){
        // Subroutine to query kmers via Jellyfish database to refine the coordinate space. 
        // Rules:
        // 1. traverse nonoverlapping 21 mers 5kb from each side per chromosome
        // 2. terminate at the first run of 21 mers above 10 copies
        // 3. if refinement doesn't work, keep the original coordinates
        KmerRepeatClassifier workhorse = new KmerRepeatClassifier(JellyfishDb);
        FastaSequenceIndex index = new FastaSequenceIndex(new File(this.origin.getPath().toAbsolutePath() + ".fai"));
        IndexedFastaSequenceFile reader = new IndexedFastaSequenceFile(this.origin.getPath(), index);
        
        
        this.mergedPlan.entrySet().stream()
                .filter(s -> this.mergedPlan.get(s.getKey()).size() > 1)
                .forEach((s) -> {
                    List<BedFastaPlan> beds = s.getValue();
                    
                    for(int x = 0; x < beds.size(); x++){
                        String contig = beds.get(x).Chr();
                        long oBegin = beds.get(x).Start();
                        long oEnd = beds.get(x).End();
                        
                        long actualE = this.origin.getChrLen(contig);
                        if(oEnd > actualE){
                            log.log(Level.WARNING, "[WARN] Actual subsequence size of chr: " + contig + " was: " + actualE);
                            oEnd = actualE;
                            beds.get(x).End();
                        }
                        
                        if(x != beds.size() -1){
                            // Work on the end coordinates if not at the end of the chromosome plan
                            long nStart = oEnd - 5000;
                            if(nStart <= 100){
                                // terminate comparison with small fragments
                                continue;
                            }                           
                            
                            int newEnd = workhorse.RefineEndCoord(reader.getSubsequenceAt(contig, nStart, oEnd)
                                    .getBaseString());
                            if(newEnd != -1){
                                int tval = (int)oEnd - ((newEnd + 1) * 21);
                                beds.get(x).setEnd(tval);
                                log.log(Level.INFO, "[REFINE] " + beds.get(x).toString() + "\tendRefine\tto: " + tval);
                            }else
                                log.log(Level.INFO, "[REFINE] " + beds.get(x).toString() + "\tendSkipped");                            
                        }
                        if(x != 0){
                            // Work on the start coordinates if not at the start of the chromosome plan
                            long nEnd = oBegin + 5000;
                            if(nEnd > oEnd){
                                // terminate comparisons that extend beyond the frag length
                                continue;
                            }
                            int newStart = workhorse.RefineStartCoord(reader.getSubsequenceAt(contig, oBegin, nEnd)
                                    .getBaseString());
                            
                            if(newStart != -1){
                                int tval = (int)oBegin - ((newStart + 1) * 21);
                                beds.get(x).setStart(tval);
                                log.log(Level.INFO, "[REFINE] " + beds.get(x).toString() + "\tstartRefine\tto: " + tval);
                            }else
                                log.log(Level.INFO, "[REFINE] " + beds.get(x).toString() + "\tstartSkipped");
                        }
                    }
                });
    }
    
    public void printOrderedListToAGP(String outfile){
        try(BufferedWriter output = Files.newBufferedWriter(Paths.get(outfile), Charset.defaultCharset())){
            // Write the modified plans
            this.mergedPlan.entrySet().stream()
                    .sorted( (s, s1) ->{
                            if(utils.MergerUtils.isNumeric(s.getKey()) 
                                    && utils.MergerUtils.isNumeric(s1.getKey())){
                                return Integer.compare(Integer.parseInt(s.getKey()), Integer.parseInt(s1.getKey()));
                            }
                            return s.getKey().compareTo(s1.getKey());})
                    .forEachOrdered((s)->{
                        List<BedFastaPlan> beds = s.getValue();
                        int prevend = 1;
                        for(int x = 0; x < beds.size() -2; x++){
                            Map<Integer, String> agp = beds.get(x).toAGP(prevend, x + 1, (beds.size() > 1));
                            prevend = agp.keySet().iterator().next();
                            String ovalue = agp.get(prevend);
                            try {
                                output.write(ovalue);
                            } catch (IOException ex) {
                                log.log(Level.SEVERE, "Error writing AGP entry to output file: " + outfile + "!", ex);
                            }
                        }
                        if(beds.size() > 1){
                            // write the last entry for the chromosome
                            Map<Integer, String> agp = beds.get(beds.size() -1).toAGP(prevend, beds.size(), false);
                            prevend = agp.keySet().iterator().next();
                            String ovalue = agp.get(prevend);
                            try {
                                output.write(ovalue);
                            } catch (IOException ex) {
                                log.log(Level.SEVERE, "Error writing AGP entry to output file: " + outfile + "!", ex);
                            }
                        }                            
                    });
            
            // Write unmodified contigs
            this.unmodified.forEach((s) ->{
                long length = this.origin.getChrLen(s);
                StringBuilder out = new StringBuilder();
                out.append(s).append("\t").append(1).append("\t").append(length).append("\t")
                        .append(1).append("\t").append("D\t").append(s).append("\t").append(1)
                        .append("\t").append(length).append("\t").append("+").append(System.lineSeparator());
                try {
                    output.write(out.toString());
                } catch (IOException ex) {
                    log.log(Level.SEVERE, "Error writing AGP entry to output file: " + outfile + "!", ex);
                }
            });
        }catch(IOException ex){
            log.log(Level.SEVERE, "Error writing AGP entry to output file: " + outfile + "!", ex);
        }
    }
    
    public void printOrderedListToFasta(String outfile){
        try(BufferedWriter output = Files.newBufferedWriter(Paths.get(outfile), Charset.defaultCharset())){
            // Trying HTSJDK subsectioning test
            FastaSequenceFile reader = new FastaSequenceFile(this.origin.getPath(), true);
            
        }catch(IOException ex){
            log.log(Level.SEVERE, "Error writing fasta entry to output file: " + outfile + "!", ex);
        }
    }
    
    private BedFastaPlan convertCoordToBed(markerCoords coord){
        BedFastaPlan value; int start = 1; int end = 1;
        if(coord.prevMarker == null || !coord.prevMarker.nChr.equals(coord.prevMarker.nChr) ){
            // Start of the chromosome
            start = 1;
        }else{
            // Average of the two coords  -1 to facilitate overlap
            start = ((coord.prevMarker.nPos + coord.nPos) / 2) - 1;
            //start = coord.prevMarker.nPos - 1;
        }
        if(coord.nextMarker == null || !coord.nextMarker.nChr.equals(coord.nextMarker.nChr)){
            // end of the chromosome
            end = (int)this.origin.getChrLen(coord.nChr);
        }else{
            // Average of the two coords + 1 to facilitate overlap
            end = ((coord.nextMarker.nPos + coord.nPos) / 2) + 1;
            //end = coord.nextMarker.nPos + 1;
        }
        value = new BedFastaPlan(coord.nChr, start, end, coord.oChr, coord.oPos);
        return value;
    }
    
    
    protected class OrgChrSort implements Comparator<markerCoords>, Serializable{

        @Override
        public int compare(markerCoords t, markerCoords t1) {
            if(utils.MergerUtils.isNumeric(t.getOChr()) 
                    && utils.MergerUtils.isNumeric(t1.getOChr())){
                return Integer.compare(Integer.parseInt(t.getOChr()), Integer.parseInt(t1.getOChr()));
            }else if(utils.MergerUtils.isNumeric(t.getOChr()) 
                    && !utils.MergerUtils.isNumeric(t1.getOChr()))
                return -1;
            else if(!utils.MergerUtils.isNumeric(t.getOChr())
                    && utils.MergerUtils.isNumeric(t1.getOChr()))
                return 1;
            return t.getOChr().compareTo(t1.getOChr());
        }
        
    }
    
    // Numerical and String sorting of chromosome names
    protected class NewChrSort implements Comparator<markerCoords>, Serializable{

        @Override
        public int compare(markerCoords t, markerCoords t1) {
            if(utils.MergerUtils.isNumeric(t.getNChr()) 
                    && utils.MergerUtils.isNumeric(t1.getNChr())){
                return Integer.compare(Integer.parseInt(t.getNChr()), Integer.parseInt(t1.getNChr()));
            }else if(utils.MergerUtils.isNumeric(t.getNChr()) 
                    && !utils.MergerUtils.isNumeric(t1.getNChr()))
                return -1;
            else if(!utils.MergerUtils.isNumeric(t.getNChr())
                    && utils.MergerUtils.isNumeric(t1.getNChr()))
                return 1;
            return t.getNChr().compareTo(t1.getNChr());
        }
    
    }
    
    protected class markerCoords{
        public final String oChr;
        public final String nChr;
        public final int oPos;
        public final int nPos;
        public int order;
        // linked list of previous and current markers to estimate breakpoints of new position coordinates
        public markerCoords prevMarker = null;
        public markerCoords nextMarker = null;
        
        public markerCoords(String line){
            line = line.trim();
            String[] segs = line.split("\t");
            String[] oSegs = ProcessReadName(segs[0]);
            
            this.oChr = oSegs[1];
            this.oPos = Integer.parseInt(oSegs[2]);
            
            this.nChr = segs[2];
            this.nPos = Integer.parseInt(segs[3]);
        }
        
        private String[] ProcessReadName(String name){
            String[] elements = name.split("\\.");
            return elements;
        }
        
        public void setPrevMarker(markerCoords prev){
            this.prevMarker = prev;
        }
        
        public void setNextMarker(markerCoords next){
            this.nextMarker = next;
        }
        
        public String getNChr(){
            return this.nChr;
        }
        
        public int getNPos(){
            return this.nPos;
        }
        
        public String getOChr(){
            return this.oChr;
        }
        
        public int getOPos(){
            return this.oPos;
        }
        
        public boolean isUnmapped(){
            return this.nChr.startsWith("*");
        }
        
        public void setOrder(int i){
            this.order = i;
        }
    }
}
