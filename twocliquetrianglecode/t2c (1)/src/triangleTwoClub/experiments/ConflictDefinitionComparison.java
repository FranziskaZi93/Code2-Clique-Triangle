package triangleTwoClub.experiments;

//import javafx.util.Pair;
import robustTwoClub.graph.RollbackRtcGraph;
import robustTwoClub.graph.RtcGraph;
import robustTwoClub.graph.Triangle;
import triangleTwoClub.ConditionChecker;
import triangleTwoClub.ITriangleTwoClubAlgorithm;
import triangleTwoClub.Output;
import triangleTwoClub.graph.ConflictGraph;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class ConflictDefinitionComparison {

    public static void main(String[] args) {
        compareConflictGraphBuildTime(200);
    }

    private static void compareConflictGraphBuildTime(int sizeCap){
        List<GraphSupplier.Pair<String, Integer>> graphs = GraphSupplier.getBachelorThesisGraphs();
        for(GraphSupplier.Pair<String, Integer> info : graphs) {
            String path = info.getKey();
            RollbackRtcGraph graph = new RollbackRtcGraph(path, info.getValue());

            if(graph.size() > sizeCap)
                continue;

            System.out.println("################################");
            System.out.println(LowerBoundComparison.getName(path) + " " + " " + graph.size() + " " +
                    graph.getEdgeCount());

            long time;
            ConflictGraph con;
            ITriangleTwoClubAlgorithm algorithm = createAlgorithm(graph, 1);

            // V1
            con = new ConflictGraph(graph);
            con.setAlgorithm(algorithm);
            con.useBetterCompatibility(false);
            time = System.currentTimeMillis();
            con.buildConflictGraph();
            time = System.currentTimeMillis() - time;
            System.out.println("V1: " + time + "ms");

            // V2
            con = new ConflictGraph(graph);
            con.setAlgorithm(algorithm);
            con.useBetterCompatibility(true);
            time = System.currentTimeMillis();
            con.buildConflictGraph();
            time = System.currentTimeMillis() - time;
            System.out.println("V2: " + time + "ms");
        }
    }

    private static void compareConflictDefinitions(){
        // TODO
    }

    private static ITriangleTwoClubAlgorithm createAlgorithm(RtcGraph g, int l){
        HashSet<Integer> marked = new HashSet<>();
        HashMap<Integer, HashSet<Triangle>> triangles = Triangle.getTriangleHashMap(g);
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
                return marked;
            }

            @Override
            public HashMap<Integer, HashSet<Triangle>> getTriangles() {
                return triangles;
            }

            @Override
            public ConditionChecker getConditionChecker() {
                return new ConditionChecker(this);
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
            public boolean isEdgeConditionEnabled() {
                return false;
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
}
