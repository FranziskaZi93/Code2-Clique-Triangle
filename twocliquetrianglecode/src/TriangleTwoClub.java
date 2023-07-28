import robustTwoClub.algorithms.TriangleTwoClubAlgorithm;
import robustTwoClub.graph.RtcGraph;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.StringTokenizer;

public class TriangleTwoClub {

    public static void main(String[] args) {

        if(args.length < 2){
            System.out.println("Triangle Two Club Algorithm Usage:");
            System.out.println("TriangleTwoClub.jar [Graph Input File] [Graph Type Number] [optional arguments]");
            System.out.println("Use -h to show help");
            System.exit(0);
        }

        if(args[0].equals("-h")){
            System.out.println("\nHelp");
            System.out.println("\nUsage:");
            System.out.println("\tTriangleTwoClub.jar [Graph Input File] [Graph Type Number] [optional arguments]");
            System.out.println("\nGraph Types:");
            System.out.println("\t 0 = DIMACS (with header line)");
            System.out.println("\t 1 = METIS");
            System.out.println("\t 2 = DIMACS clq");
            System.out.println("\t 3 = edges");
            System.out.println("\nBranching Methods:");
            System.out.println(" 0 = Single vertex branching (default)");
            System.out.println(" 1 = Incompatible vertices branching");
            System.out.println("\n Arguments:");
            System.out.println("-h \tshow this dialog");
            System.out.println("-o \tprint result vertex names in console");
            System.out.println("-s \tsilent run (no algorithm output)");
            System.out.println("-t [seconds] \ttime limit in seconds as float (5.0 = 5 seconds)");
            System.out.println("-f [path] \tset output file");
            System.out.println("-b [method] \tsets the branching method");
            System.out.println("-bs [steps] \tuse reduction every [steps] branching depth (valid: 1-10)");
            System.out.println("-noVC \tdisables vertex cover rule");
            System.out.println("-vcs [steps] \tuse vertex cover rule every [steps] branching depth (valid: 1-10)");
            System.out.println("-noNCR \tdisables no choice rule");
            System.exit(0);
        }

        if(!new File(args[0]).exists()){
            System.out.println("Input file not found.");
            System.exit(1);
        }
        int graphType = 0;

        try{
            graphType = Integer.parseInt(args[1]);
            assert graphType >= 0 && graphType < 3;
        }catch(Exception e){
            System.out.println("Invalid graph type: " + args[1]);
            System.exit(1);
        }

        // determine settings
        String graphName = new File(args[0]).getName();
        boolean printResultVertices = false;
        boolean silentRun = false;
        boolean useVC = true;
        boolean useNCR = true;
        float timeLimit = 0.0f;
        String outputPath = null;
        TriangleTwoClubAlgorithm.Branching method = TriangleTwoClubAlgorithm.Branching.SINGLE_VERTEX;
        int branchingReductionSteps = 1;
        int vertexCoverRuleSteps = 1;

        int innerI = 0;
        try {
            for (int i = 0; i < args.length; i++) {
                innerI = i;
                if (args[i].equals("-o")) printResultVertices = true;
                if (args[i].equals("-s")) silentRun = true;
                if (args[i].toLowerCase().equals("-novc")) useVC = false;
                if (args[i].toLowerCase().equals("-noncr")) useNCR = false;
                if (args[i].equals("-t")) {
                    timeLimit = Float.parseFloat(args[++i]);
                }
                if (args[i].equals("-bs")) {
                    branchingReductionSteps = Integer.parseInt(args[++i]);
                }
                if (args[i].equals("-vcs")) {
                    vertexCoverRuleSteps = Integer.parseInt(args[++i]);
                }
                if (args[i].equals("-f")) {
                    outputPath = args[++i];
                    File f = new File(outputPath);
                    if (f.exists() && !f.canWrite()) {
                        System.out.println("Can not write to file \"" + outputPath + "\"");
                        System.exit(0);
                    }
                }
                if (args[i].equals("-b")) {
                    if (args[++i].equals("0"))
                        method = TriangleTwoClubAlgorithm.Branching.SINGLE_VERTEX;
                    else if (args[i].equals("1"))
                        method = TriangleTwoClubAlgorithm.Branching.INCOMPATIBLE_VERTICES;
                    else {
                        System.out.println("Invalid branching method: " + args[i]);
                        System.exit(0);
                    }
                }
            }
        }catch (ArrayIndexOutOfBoundsException e){
            System.out.println("Parameter " + args[innerI] + " requires an additional argument");
            System.exit(0);
        }catch (NumberFormatException e){
            System.out.println(args[innerI] + " is not a valid time");
            System.exit(0);
        }



        // run algorithm
        RtcGraph graph = buildGraphFromFile(args[0], graphType);
        TriangleTwoClubAlgorithm.setBranchReductionSteps(branchingReductionSteps);
        TriangleTwoClubAlgorithm.setTimeLimit(timeLimit);
        TriangleTwoClubAlgorithm.setBranchingMethod(method);
        TriangleTwoClubAlgorithm.useVertexCoverRule(useVC);
        TriangleTwoClubAlgorithm.useNoChoiceRule(useNCR);
        TriangleTwoClubAlgorithm.setVertexCoverRuleSteps(vertexCoverRuleSteps);
        if(outputPath != null) TriangleTwoClubAlgorithm.setOutputFile(outputPath);
        RtcGraph result = TriangleTwoClubAlgorithm.run(graph, graphName, !silentRun);

        // output
        if(printResultVertices){
            List<String> names = new LinkedList<>();
            for(int v : result.getVertices()) names.add(result.getVertexName(v));
            System.out.println(names);
        }

    }

    public static RtcGraph buildGraphFromFile(String filename, int type){
        RtcGraph graph = new RtcGraph();
        HashMap<String, Integer> seenIds = new HashMap<String, Integer>();	// only required during construction
        File file = new File(filename);

        try {
            final LineNumberReader reader = new LineNumberReader(new FileReader(file));
            String str;
            int id = 0;
            if (type == 0) {		// Parse DIMACS file
                reader.readLine();  //read first line containing the number of vertices and edges
                while ((str = reader.readLine()) != null)
                {
                    str = str.trim();										// trim away whitespace at either end of line
                    if (!str.startsWith("#") && !str.startsWith("%")) {								// skip comment lines
                        StringTokenizer tokens = new StringTokenizer(str);
                        if (tokens != null && tokens.countTokens() > 1) {	// only consider well-formed lines
                            String vertexA = tokens.nextToken();
                            String vertexB = tokens.nextToken();
                            if (!seenIds.containsKey(vertexA)) {		// add vertex 0 if never seen before
                                seenIds.put(vertexA, id);
                                graph.addVertex(id, vertexA);
                                id++;
                            }
                            if (!seenIds.containsKey(vertexB)) {		// add vertex 1 if never seen before
                                seenIds.put(vertexB, id);
                                graph.addVertex(id, vertexB);
                                id++;
                            }
                            // add edge to both adjacency lists
                            graph.addEdge(seenIds.get(vertexA), seenIds.get(vertexB));
                        }
                    }
                }
            } else if (type == 1) {		// parse METIS file
                StringTokenizer header = new StringTokenizer(reader.readLine());
                int vertices = Integer.parseInt(header.nextToken());
                int edges = Integer.parseInt(header.nextToken());
                boolean weights = header.hasMoreTokens() && !header.nextToken().equals("0");
                for (int v = 0; v < vertices; v++) {
                    graph.addVertex(v, Integer.toString(v+1));
                }
                int edgesCounted = 0;
                int currentVertex = 0;
                while ((str = reader.readLine()) != null)
                {
                    if (!str.startsWith("#") && !str.startsWith("%")) {				// skip comment lines
                        StringTokenizer tokens = new StringTokenizer(str);
                        while (tokens.hasMoreTokens()) {
                            int v = Integer.parseInt(tokens.nextToken()) - 1;
                            if (!graph.adjacent(currentVertex, v)) {
                                graph.addEdge(currentVertex, v);
                                edgesCounted++;
                            }
                            if (weights) tokens.nextToken(); // Throw away weight token
                        }
                        currentVertex++;
                    }
                }
                if (edges != edgesCounted) {
                    System.out.println("WARNING: Wrong edge count while loading METIS file "+filename+".");
                    System.out.println("Edges stated: "+edges+"   Edges counted: "+edgesCounted+"\n");
                }
            } else if (type == 2) {
                while ((str = reader.readLine()) != null)
                {
                    StringTokenizer tokens = new StringTokenizer(str);
                    String lineType = tokens.nextToken();
                    if (lineType.equals("c")) continue;
                    if (lineType.equals("p")) {
                        tokens.nextToken();
                        int vertices = Integer.parseInt(tokens.nextToken());
                        for (int v = 0; v < vertices; v++) {
                            graph.addVertex(v, Integer.toString(v+1));
                        }
                    }
                    if (lineType.equals("e")) {
                        graph.addEdge(Integer.parseInt(tokens.nextToken()) - 1, Integer.parseInt(tokens.nextToken()) - 1);
                    }
                }
            } else if (type == 3) {		// Edges only
                while ((str = reader.readLine()) != null)
                {
                    str = str.trim();										// trim away whitespace at either end of line
                    if (!str.startsWith("#") && !str.startsWith("%")) {								// skip comment lines
                        StringTokenizer tokens = new StringTokenizer(str);
                        if (tokens != null && tokens.countTokens() > 1) {	// only consider well-formed lines
                            String vertexA = tokens.nextToken();
                            String vertexB = tokens.nextToken();
                            if (!seenIds.containsKey(vertexA)) {		// add vertex 0 if never seen before
                                seenIds.put(vertexA, id);
                                graph.addVertex(id, vertexA);
                                id++;
                            }
                            if (!seenIds.containsKey(vertexB)) {		// add vertex 1 if never seen before
                                seenIds.put(vertexB, id);
                                graph.addVertex(id, vertexB);
                                id++;
                            }
                            // add edge to both adjacency lists
                            graph.addEdge(seenIds.get(vertexA), seenIds.get(vertexB));
                        }
                    }
                }
            }
            reader.close();
        } catch (IOException e) {
            System.out.println("Could not locate input file '"+filename+"'.");
            System.exit(0);
        }
        return graph;
    }

}
