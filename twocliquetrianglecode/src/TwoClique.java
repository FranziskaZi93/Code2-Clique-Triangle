import robustTwoClub.graph.RollbackRtcGraph;
import robustTwoClub.graph.RtcGraph;
import triangleTwoClub.Algorithm;
import triangleTwoClub.LowerBound;
import triangleTwoClub.Output;
import triangleTwoClub.graph.T2CGraph;

import java.util.Scanner;

public class TwoClique {
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);  // Create a Scanner object
        System.out.println("Add your comand:");

        String strComand = scanner.nextLine();  // Read user input
        String[] strArray = null;
        strArray = strComand.split(" ");

        runTriangleTwoClub(strArray);

    }


    private static void runTriangleTwoClub(String[] args){
        // args[1] = graph file
        // args[2] = graph type
        if(args.length < 3){
            return;
        }


        int graphType;
        try {
            graphType = Integer.parseInt(args[2]);
            if(graphType < 0 || graphType > 4)
                throw new RuntimeException();
        }catch(Exception e){
            System.out.println("Invalid graph type: " + args[2]);
            return;
        }

        int ell = 1, timeLimit = -1;
        String output = null;
        boolean edgeMode = false, silentMode = false, singleLine = false, addSolution = false;
        LowerBound.LowerBoundVariant lbVariant = LowerBound.LowerBoundVariant.LB1;
        T2CGraph.ConflictGraphType conflictGraphType = T2CGraph.ConflictGraphType.SIMPLE;
        for(int i = 3; i < args.length; i++){
            if(args[i].equals("-o"))
                output = args[++i];

            if(args[i].equals("-t"))
                timeLimit = Integer.parseInt(args[++i]);

            if(args[i].equals("-l"))
                ell = Integer.parseInt(args[++i]);

            if(args[i].equals("-lb"))
                switch(args[++i]){
                    case "b": lbVariant = LowerBound.LowerBoundVariant.BASIC; break;
                    case "0": lbVariant = LowerBound.LowerBoundVariant.DISABLED; break;
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

    private static String getName(String path){
        if(path.contains("\\"))
            return path.substring(path.lastIndexOf("\\") + 1);
        else
            return path.substring(path.lastIndexOf("/") + 1);
    }
}
