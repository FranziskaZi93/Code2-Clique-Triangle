package triangleTwoClub;

import robustTwoClub.graph.RtcGraph;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Locale;

public class Output {

    /**
     * Specifies the format of stored values.
     */
    public enum DataFormat{PLAIN, MILLISECONDS, NANOSECONDS, PERCENTAGE, INTEGER}
    public enum OutputFormat{LINE_SEPARATED, SYMBOL_SEPARATED}

    public static long consolePrintDelay = 10000L; // minimum delay between console prints in milliseconds
    private static char separatorSymbol = ';';
    private final DecimalFormat timeFormat;
    private final DecimalFormat decimalFormat;

    private OutputFormat outputFormat;
    private LinkedList<String> information;
    private HashMap<String, Double> statistics;
    private HashMap<String, DataFormat> formats;
    private long lastConsolePrint;

    public Output(){
        reset();
        outputFormat = OutputFormat.LINE_SEPARATED;
        DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.getDefault());
        symbols.setDecimalSeparator('.');
        symbols.setGroupingSeparator(',');
        timeFormat = new DecimalFormat("0.000");	// Format for time output;
        timeFormat.setDecimalFormatSymbols(symbols);
        decimalFormat = new DecimalFormat("00.00"); // Format for progress output
    }

    /**
     * Adds an information for the output. An information (e.g. the name of the graph) gets printed in front of the
     * statistics in the output file.
     * For example: key = Name, value = TestGraph produces "Name=TestGraph"
     * Adding a value for the same key twice overrides the first one.
     * @param key the key of the information
     * @param value the value of the information
     */
    public void addInformation(String key, String value){
        String info = key + "=" + value.trim();

        String foundEntry = null;
        for(String s : information)
            if(s.startsWith(key + "=")){
                foundEntry = s;
                break;
            }
        if(foundEntry != null)
            information.remove(foundEntry);

        information.add(info);
    }

    /**
     * Logs a value. If this key has been used before, then the new value is added to the current one.
     * @param key the name of the logged value
     * @param value the value to log
     */
    public void log(String key, double value){
        log(key, value, true);
    }

    /**
     * Logs a value. If this key has been used before, then the new value is added to the current one.
     * @param key the name of the logged value
     * @param value the value to log
     */
    public void log(String key, double value, boolean accumulate){
        if(statistics.containsKey(key) && accumulate)
            statistics.put(key, statistics.get(key) + value);
        else
            statistics.put(key, value);
    }

    /**
     * Specifies the format of logged values.
     * @param key the key of a value
     * @param format the output format of the key
     */
    public void setFormat(String key, DataFormat format){
        formats.put(key, format);
    }

    /**
     * Specifies the format of the output file.
     * @param f the format
     */
    public void setOutputFormat(OutputFormat f){
        this.outputFormat = f;
    }

    /**
     * Prints the current status of the algorithms run to the console.
     */
    public void printStatusToConsole(){
        long time = System.currentTimeMillis();
        if(time < lastConsolePrint + consolePrintDelay)
            return;
        try{
            System.out.print("start point " + statistics.get("CurrentStartPoint").intValue());
            System.out.print("/" + statistics.get("StartPointCount").intValue());
            System.out.print(" - " + statistics.get("BranchCount").intValue() + " branches - ");
            System.out.println(" max depth " + statistics.get("MaxBranchDepth").intValue());
        }catch(Exception e){
            // a failed status print shall not stop the algorithms execution
        }
        lastConsolePrint = time;
    }

    /**
     * Writes all statistics to a file.
     * @param path the path of the output file.
     */
    public void writeOutputFile(String path){
        File f = new File(path);
        try {
            PrintWriter writer = new PrintWriter(new FileWriter(f, true), true);

            if(outputFormat == OutputFormat.LINE_SEPARATED){
                // print information first
                for(String s : information)
                    writer.println(s);

                // print statistics
                for(String key : statistics.keySet()){
                    String output = createOutputLine(key).trim();
                    if(!output.isEmpty())
                        writer.println(output);
                }
            }
            else if(outputFormat == OutputFormat.SYMBOL_SEPARATED){
                for(String s : information){
                    writer.print(s);
                    writer.print(separatorSymbol);
                }

                for(String key : statistics.keySet()){
                    String output = createOutputLine(key).trim();
                    if(!output.isEmpty()){
                        writer.print(output);
                        writer.print(separatorSymbol);
                    }
                }
                writer.println();
            }

            writer.flush();
            writer.close();
        } catch (IOException e) {
            System.out.println("Failed to write output to file \"" + path + "\".");
            e.printStackTrace();
        }
    }

    /**
     * Deletes all stored information, logged statistics and formats.
     */
    public void reset(){
        information = new LinkedList<>();
        statistics = new HashMap<>();
        formats = new HashMap<>();
        lastConsolePrint = System.currentTimeMillis();
    }

    /**
     * Get a logged statistic.
     * @param key the key of the logged statistic
     * @return the value associated with the key
     */
    public Double getLoggedValue(String key){
        return statistics.get(key);
    }

    /**
     * Get the value of a stored information.
     * @param key the key of the information
     * @return the value assigned to the key, an empty string if the key does not exist
     */
    public String getInformation(String key){
        for(String s : information)
            if(s.startsWith(key + "="))
                return s.substring(key.length() + 1);
         return "";
    }

    /**
     * Logs the vertices of a graph in the output.
     * @param key the key for the information
     * @param graph a graph object
     */
    public void logGraph(String key, RtcGraph graph){
        StringBuilder graphString = new StringBuilder();
        graphString.append("{");

        if(graph != null) {
            int n = graph.size(), i = 1;
            for (Integer v : graph.getVertices()) {
                graphString.append(graph.getVertexName(v));
                if (i < n)
                    graphString.append(',');
                i++;
            }
        }

        graphString.append("}");
        addInformation(key, graphString.toString());
    }

    private String createOutputLine(String key){
        if(!statistics.containsKey(key))
            return key.trim() + "=???";

        double value = statistics.get(key);
        DataFormat format = formats.get(key);

        if(format == null || format == DataFormat.PLAIN)
            return key.trim() + "=" + value;

        if(format == DataFormat.INTEGER)
            return key.trim() + "=" + ((int) value);

        if(format == DataFormat.PERCENTAGE)
            return key.trim() + "=" + decimalFormat.format(value) + "%";

        if(format == DataFormat.MILLISECONDS)
            return key.trim() + "=" + timeFormat.format(value / 1000.0);

        // only nanoseconds is left
        return key.trim() + "=" + timeFormat.format(value / 1000000000.0);
    }
}
