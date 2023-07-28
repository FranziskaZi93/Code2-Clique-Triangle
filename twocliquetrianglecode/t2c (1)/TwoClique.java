import robustTwoClub.graph.RollbackRtcGraph;
import robustTwoClub.graph.RtcGraph;
import triangleTwoClub.Algorithm;
import triangleTwoClub.LowerBound;
import triangleTwoClub.Output;
import triangleTwoClub.graph.T2CGraph;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class TwoClique {




        public static void main(String[] args) {
            if(args.length == 0 || args[0].equals("-h")){
              //  printHelp();
                return;
            }//* -help wird aufgerufen

            if(args[0].equals("l")) //* wofür steht l vergleich lower bounds
                runLowerBoundComparison(args);
            else if(args[0].equals("t")) // wofür steht triangle two clubs aufgerufen
                runTriangleTwoClub(args);
        }

        private static void runTriangleTwoClub(String[] args){
            // args[1] = graph file
            // args[2] = graph type
            if(args.length < 3){ //* das bedeutet dann kann es kein Dreieck werden
               // printTriangleHelp();
                return;
            }


            int graphType;
            try {
                graphType = Integer.parseInt(args[2]);
                if(graphType < 0 || graphType > 4)
                    throw new RuntimeException();  //* warum ? funktionen graphtypen, welche Graphtypen wo sind sie,
            }catch(Exception e){
                System.out.println("Invalid graph type: " + args[2]);
              //  printGraphTypes();
                return;
            }

            int ell = 1,  timeLimit = -1; //*kein zeitlimit, wird nicht benutzt,
            String output = null;
            boolean edgeMode = false, silentMode = false, singleLine = false, addSolution = false; //*was heißt das alles, warum wird edgeMode auf False gesetzt, was bedeutet silentMode; was bedeutet singleLine, warum ist addSolution falsch
            LowerBound.LowerBoundVariant lbVariant = LowerBound.LowerBoundVariant.LB1;
            T2CGraph.ConflictGraphType conflictGraphType = T2CGraph.ConflictGraphType.SIMPLE;
            for(int i = 3; i < args.length; i++){
                if(args[i].equals("-o")) //* wofür steht o
                    output = args[++i];

                if(args[i].equals("-t"))
                    timeLimit = Integer.parseInt(args[++i]);

                if(args[i].equals("-l"))
                    ell = Integer.parseInt(args[++i]);

                if(args[i].equals("-lb")) //paramter lower bound
                    switch(args[++i]){
                        case "b": lbVariant = LowerBound.LowerBoundVariant.BASIC; break;
                        case "0": lbVariant = LowerBound.LowerBoundVariant.DISABLED; break; //*das hier setzen
                        case "1": lbVariant = LowerBound.LowerBoundVariant.LB1; break;
                        case "2": lbVariant = LowerBound.LowerBoundVariant.LB2; break;
                    }

                if(args[i].equals("-cg"))
                    switch(args[++i]){
                        case "0": conflictGraphType = T2CGraph.ConflictGraphType.DISABLED; break;
                        case "1": conflictGraphType = T2CGraph.ConflictGraphType.SIMPLE; break;
                        case "2": conflictGraphType = T2CGraph.ConflictGraphType.ADVANCED; break;
                    }



                if(args[i].equals("-e"))
                    edgeMode = true;

                if(args[i].equals("-s"))
                    silentMode = true;

                if(args[i].equals("-SL"))
                    singleLine = true;

                if(args[i].equals("-og"))
                    addSolution = true;

                //* was ist das alles
            }

            Algorithm t2c = new Algorithm();
            t2c.setL(ell);
            t2c.enableEdgeCondition(edgeMode);
            t2c.setTimeLimit(timeLimit * 1000L);
            t2c.setLowerBoundVariant(lbVariant);
            t2c.setConflictGraphType(conflictGraphType);
            t2c.addSolutionToOutput(addSolution);
            t2c.getOutput().setOutputFormat(
                    singleLine ? Output.OutputFormat.SYMBOL_SEPARATED : Output.OutputFormat.LINE_SEPARATED
            );

            RollbackRtcGraph graph = new RollbackRtcGraph(args[1], graphType);
            RtcGraph solution = t2c.run(graph, getName(args[1]), output, !silentMode);
            // TODO: make it possible to save solution as graph
        }

        private static void runLowerBoundComparison(String[] args){
            // args[1] = graph file
            // args[2] = graph type
            if(args.length < 3){
             //   printHelpLowerBoundComparison();
                return;
            }


            int graphType;
            try {
                graphType = Integer.parseInt(args[2]);
                if(graphType < 0 || graphType > 4)
                    throw new RuntimeException();
            }catch(Exception e){
                System.out.println("Invalid graph type: " + args[2]);
               // printGraphTypes();
                return;
            }


            int cap = 100, timeLimit = -1;
            String output = null;
            boolean edgeMode = false;
            for(int i = 3; i < args.length; i++){
                if(args[i].equals("-o"))
                    output = args[++i];

                if(args[i].equals("-t"))
                    timeLimit = Integer.parseInt(args[++i]);

                if(args[i].equals("-c"))
                    cap = Integer.parseInt(args[++i]);

                if(args[i].equals("-e"))
                    edgeMode = true;
            }
/**
            StringBuilder results = LowerBoundComparison.runComparison(args[1], graphType, cap, timeLimit, edgeMode);
            if(output == null){
                System.out.println(results.toString());
            }else{
                PrintWriter writer;
                try {
                    FileWriter fw = new FileWriter(new File(output));
                    writer = new PrintWriter(fw, true);
                    writer.print(results.toString());
                    writer.flush();
                    writer.close();
                } catch (IOException e) {
                    e.printStackTrace();
                    return;
                }
            }
 */
        }

 /**

        private static void printHelp(){
            System.out.println("\n2-Club Algorithms");
            System.out.println("Usage: TwoClubs.jar [type]");
            System.out.println("Types:");
            System.out.println("\t t - triangle-two-clubs");
            System.out.println("\t l - triangle-two-club lower bound comparison");
            // TODO
        }
*/
 /**
        private static void printTriangleHelp(){
            System.out.println("Triangle Two Club Usage:");
            System.out.println("TwoClubs.jar t [input graph] [graph type] [arguments]\n");
            System.out.println("Arguments:");
            System.out.println("\t -o [path] - print result information to file");
            System.out.println("\t -og - adds the vertices of the solution to the result information file");
            System.out.println("\t -t [limit] - set a time limit in seconds");
            System.out.println("\t -l [value] - set a value for L (default is 1)");
            System.out.println("\t -lb [type] - choose a lower bound variant (default=1)");
            System.out.println("\t\t -> B = basic, only works for L=1");
            System.out.println("\t\t -> 0 = lower bound disabled");
            System.out.println("\t\t -> 1 = works for all L (LB1)");
            System.out.println("\t\t -> 2 = works for all L, deeper search (LB2)");
            System.out.println("\t -cg [type] - choose the conflict graph variant (default=1)");
            System.out.println("\t\t -> 0 = conflict graph disabled");
            System.out.println("\t\t -> 1 = simple (distance only conflict)");
            System.out.println("\t\t -> 2 = advanced conflict definition");
            System.out.println("\t -e - enable edge-mode");
            System.out.println("\t -s - enable silent mode (no console prints)");
            System.out.println("\t -SL - single line output in output file (separated by symbol)");
        }

        private static void printHelpLowerBoundComparison(){
            System.out.println("Lower Bound Comparison Usage:");
            System.out.println("TwoClubs.jar l [input graph] [graph type] [arguments]\n");
            System.out.println("Arguments:");
            System.out.println("\t -o [path] - print result to file");
            System.out.println("\t -t [limit] - set a time limit in milliseconds");
            System.out.println("\t -c [cap] - the maximum value for L");
            System.out.println("\t -e - enable edge-mode");
        }

        private static void printGraphTypes(){
            System.out.println("Valid Graph Types:");
            System.out.println("\t 0 = .dimacs, .edges, .mtx");
            System.out.println("\t 1 = .graph");
            System.out.println("\t 2 = .clq, .col");
        }
*/
        private static String getName(String path){
            if(path.contains("\\"))
                return path.substring(path.lastIndexOf("\\") + 1);
            else
                return path.substring(path.lastIndexOf("/") + 1);
        }

    }


