package triangleTwoClub;

import robustTwoClub.graph.RtcGraph;

public class UpperBound {

    private ITriangleTwoClubAlgorithm algorithm;
    private Output output;
    private static String upperBoundCallCount = "UpperBoundCalls";
    private static String timeMeasurementName = "UpperBoundTime[s]";

    public UpperBound(ITriangleTwoClubAlgorithm algorithm){
        this.algorithm = algorithm;
        output = algorithm.getOutput();
        output.setFormat(timeMeasurementName, Output.DataFormat.MILLISECONDS);
        output.setFormat(upperBoundCallCount, Output.DataFormat.INTEGER);
    }

    /**
     * Computed the upper bound for a vertex.
     * @param v the vertex id
     * @return the upper bound for vertex v
     */
    public int getUpperBound(int v){
        output.log(upperBoundCallCount, 1);
        long time = System.currentTimeMillis();
        int u = calculateUpperBound(algorithm.getGraph(), v);
        time = System.currentTimeMillis() - time;
        output.log(timeMeasurementName, time);
        return u;
    }

    /**
     * Computed the upper bound for a vertex in a graph.
     * @param graph the graph
     * @param v the vertex id
     * @return the upper bound for vertex v
     */
    public static int calculateUpperBound(RtcGraph graph, int v){
        int sum = 0;
        for(Integer w : graph.getNeighbors(v)) sum += graph.getNeighbors(w).size();
        return sum + 1;
    }

}
