package triangleTwoClub.experiments;

import robustTwoClub.graph.RollbackRtcGraph;
import robustTwoClub.graph.RtcGraph;
import robustTwoClub.graph.Triangle;
import triangleTwoClub.ConditionChecker;
import triangleTwoClub.ITriangleTwoClubAlgorithm;
import triangleTwoClub.LowerBound;
import triangleTwoClub.Output;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class LowerBoundComparison {

    public static void main(String[] args) {
        //simpleComparison();
        //GraphSupplier.checkPaths();
        compareLB1andLB2(70, 12, 2, 600000L);
    }

    public static StringBuilder runComparison(String path, int graphType, int cap, long timeLimit, boolean edgeMode){
        StringBuilder log = new StringBuilder();
        System.out.println("starting ...");

        RollbackRtcGraph graph = new RollbackRtcGraph(path, graphType);
        HashMap<Integer, HashSet<Triangle>> triangles = Triangle.getTriangleHashMap(graph);
        log.append("# " + getName(path)+ " " + graph.size() + " " + graph.getEdgeCount() + "\n#\n");

        long time;
        boolean abortLB1 = false, abortLB2 = false;
        int zeroCountLB1 = 0, zeroCountLB2 = 0;
        for(int l = 1; l <= cap; ){
            System.out.println("L = " + l);
            log.append("L=").append(l).append("\n");
            ITriangleTwoClubAlgorithm algorithm = createAlgorithm(graph, triangles, l, edgeMode);
            LowerBound lb;
            RtcGraph solution;

            // LB1
            if(abortLB1){
                log.append("LB1 aborted\n");
            }else{
                time = System.currentTimeMillis();
                lb = new LowerBound(algorithm, LowerBound.LowerBoundVariant.LB1);
                lb.setTimeLimit(timeLimit);
                solution = lb.getLowerBoundSolution();
                time = System.currentTimeMillis() - time;
                int solutionSize = solution == null ? 0 : solution.size();
                log.append("LB1 " + solutionSize + " " + time + "ms\n");
                if(solutionSize == 0){
                    if(++zeroCountLB1 >= 5)
                        abortLB1 = true;
                }
            }

            // LB2
            if(abortLB2){
                log.append("LB2 aborted\n");
            }else{
                time = System.currentTimeMillis();
                lb = new LowerBound(algorithm, LowerBound.LowerBoundVariant.LB2);
                lb.setTimeLimit(timeLimit);
                solution = lb.getLowerBoundSolution();
                time = System.currentTimeMillis() - time;
                int solutionSize = solution == null ? 0 : solution.size();
                log.append("LB2 " + solutionSize + " " + time + "ms\n");
                if(solutionSize == 0){
                    if(++zeroCountLB2 >= 5)
                        abortLB2 = true;
                }
            }

            log.append("\n#\n");
            l++;
            if(l > 15)
                l++;
        }

        return log;
    }

    private static void simpleComparison(){
        List<GraphSupplier.Pair<String, RtcGraph>> graphs = GraphSupplier.getPaperGraphs();

        LowerBound lb;
        RtcGraph graph;
        HashMap<Integer, HashSet<Triangle>> triangles;
        long time;
        RtcGraph solution;
        int size;
        for(int l = 1; l < 5; l++) {
            System.out.println("#####################################");
            System.out.println("L=" + l);
            for (int i = 0; i < graphs.size(); i++) {
                System.out.println("Graph=" + graphs.get(i).getKey());
                graph = graphs.get(i).getValue();
                triangles = Triangle.getTriangleHashMap(graph);

                // basic
                if(l == 1){
                    lb = new LowerBound(createAlgorithm(graph, triangles, l, false), null,
                            LowerBound.LowerBoundVariant.BASIC);
                    time = System.currentTimeMillis();
                    solution = lb.getLowerBoundSolution();
                    time = System.currentTimeMillis() - time;
                    size = solution == null ? 0 : solution.size();
                    System.out.println("Basic\t" + size + "\t" + time + "ms");
                }

                // LB1
                lb = new LowerBound(createAlgorithm(graph, triangles, l, false), null,
                        LowerBound.LowerBoundVariant.LB1);
                time = System.currentTimeMillis();
                solution = lb.getLowerBoundSolution();
                time = System.currentTimeMillis() - time;
                size = solution == null ? 0 : solution.size();
                System.out.println("LB1\t\t" + size + "\t" + time + "ms");

                // LB2
                lb = new LowerBound(createAlgorithm(graph, triangles, l, false), null,
                        LowerBound.LowerBoundVariant.LB2);
                time = System.currentTimeMillis();
                solution = lb.getLowerBoundSolution();
                time = System.currentTimeMillis() - time;
                size = solution == null ? 0 : solution.size();
                System.out.println("LB2\t\t" + size + "\t" + time + "ms");

                System.out.println();
            }
        }
        System.out.println("#####################################");
    }

    private static void compareLB1andLB2(int sizeCap, int lCap, int lSteps, long timeLimit){
        String outputPath = "C:\\Users\\phili\\Documents\\Uni\\Job Triangl2Club\\Stuff\\comparison.txt";
        PrintWriter writer;
        try {
            FileWriter fw = new FileWriter(new File(outputPath));
            writer = new PrintWriter(fw, true);
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

        writer.println("# Format: LBx [size of solution] [time to compute]");
        writer.println("# time limit = " + timeLimit + "ms");
        List<GraphSupplier.Pair<String, Integer>> graphs = GraphSupplier.getBachelorThesisGraphs();
        for(GraphSupplier.Pair<String, Integer> info : graphs){
            String path = info.getKey();
            RtcGraph graph = new RtcGraph(path, info.getValue());

            if(graph.size() > sizeCap)
                continue;

            HashMap<Integer, HashSet<Triangle>> triangles = Triangle.getTriangleHashMap(graph);

            writer.println("################################");
            writer.println(getName(path) + " " + " " + graph.size() + " " + graph.getEdgeCount());

            boolean abortLB1 = false, abortLB2 = false;
            int zeroCountLB1 = 0, zeroCountLB2 = 0;
            for(int l = 1, stepCount = 0; l <= lCap; ){
                writer.println("L=" + l);
                ITriangleTwoClubAlgorithm algorithm = createAlgorithm(graph, triangles, l, true);
                RtcGraph solution;
                LowerBound lb;
                int solutionSize;
                long time;

                // LB 1
                if(abortLB1){
                    writer.println("LB1 aborted");
                }else{
                    time = System.currentTimeMillis();
                    lb = new LowerBound(algorithm, LowerBound.LowerBoundVariant.LB1);
                    lb.setTimeLimit(timeLimit);
                    solution = lb.getLowerBoundSolution();
                    time = System.currentTimeMillis() - time;
                    solutionSize = solution == null ? 0 : solution.size();
                    writer.println("LB1 " + solutionSize + " " + time + "ms");
                    if(solutionSize == 0){
                        if(++zeroCountLB1 >= 5)
                            abortLB1 = true;
                    }
                }


                // LB 2
                if(abortLB2){
                    writer.println("LB1 aborted");
                }else {
                    time = System.currentTimeMillis();
                    lb = new LowerBound(algorithm, LowerBound.LowerBoundVariant.LB2);
                    lb.setTimeLimit(timeLimit);
                    solution = lb.getLowerBoundSolution();
                    time = System.currentTimeMillis() - time;
                    solutionSize = solution == null ? 0 : solution.size();
                    writer.println("LB2 " + solutionSize + " " + time + "ms");
                    writer.println();
                    if(solutionSize == 0){
                        if(++zeroCountLB2 >= 5)
                            abortLB2 = true;
                    }
                }

                if(abortLB1 && abortLB2)
                    break;

                if(l <= 15){
                    l++;
                }else{
                    if(++stepCount >= 5){
                        stepCount = 0;
                        lSteps *= 2;
                    }else{
                        l += lSteps;
                    }
                }

            }
            System.out.println(getName(path) + " finished");
        }
    }

    private static List<HashMap<Integer, HashSet<Triangle>>> calculateTriangles(List<GraphSupplier.Pair<String,
            RtcGraph>>  graphs){
        ArrayList<HashMap<Integer, HashSet<Triangle>>> triangleMaps = new ArrayList<>(graphs.size());
        for(int i = 0; i < graphs.size(); i++){
            triangleMaps.add(i, Triangle.getTriangleHashMap(graphs.get(i).getValue()));
        }
        return triangleMaps;
    }

    private static ITriangleTwoClubAlgorithm createAlgorithm(RtcGraph g, HashMap<Integer, HashSet<Triangle>> t, int l,
                                                             boolean edges){
        return new ITriangleTwoClubAlgorithm() {
            @Override
            public RtcGraph getGraph() {
                return g;
            }

            @Override
            public RtcGraph getConflictGraph() {
                return null;
            }

            @Override
            public HashSet<Integer> getMarked() {
                return null;
            }

            @Override
            public HashMap<Integer, HashSet<Triangle>> getTriangles() {
                return t;
            }

            @Override
            public ConditionChecker getConditionChecker() {
                return null;
            }

            @Override
            public Output getOutput() {
                return null;
            }

            @Override
            public int getL() {
                return l;
            }

            @Override
            public int getBestSolutionSize() {
                return 0;
            }

            @Override
            public boolean isEdgeConditionEnabled(){
                return edges;
            }

            @Override
            public void abortBranch() {

            }

            @Override
            public void deleteVertex(int v) {

            }

            @Override
            public void deleteEdge(int v, int w) {

            }

            @Override
            public void markVertex(int v) {

            }
        };
    }

    static String getName(String path){
        if(path.contains("\\"))
            return path.substring(path.lastIndexOf("\\") + 1);
        else
            return path.substring(path.lastIndexOf("/") + 1);
    }
}