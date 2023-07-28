package triangleTwoClub.dataReduction;

import robustTwoClub.algorithms.TriangleTwoClubAlgorithm;
import triangleTwoClub.ITriangleTwoClubAlgorithm;
import triangleTwoClub.Output;

import java.util.HashSet;
import java.lang.Math;

/**
 *  A Data Reduction Rule to make the graph smaller before the actual algorithm starts.
 *  Based on the value for L, vertices with a degree smaller than a certain threshold get deleted.
 *  This rule does not require data structures provided by the algorithm.
 */
public class VertexDegreeRule extends DeleteRule {

    private int threshold;
    private String logName;

    public VertexDegreeRule(ITriangleTwoClubAlgorithm algorithm, String logName) {
        super(algorithm);
        this.logName = logName;
        this.threshold = 1;
        if(algorithm.isEdgeConditionEnabled()) {
            this.threshold = algorithm.getL() + 1;
        }else{
            double d = 0.5d + Math.sqrt(0.25d + 2 * algorithm.getL());
            this.threshold = (int) Math.ceil(d);
        }
        String thresholdKey = getRuleName() + "Threshold";
        algorithm.getOutput().setFormat(thresholdKey, Output.DataFormat.INTEGER);
        algorithm.getOutput().log(thresholdKey, this.threshold);
    }

    public VertexDegreeRule(ITriangleTwoClubAlgorithm algorithm){
        this(algorithm, "VertexDegreeRule");
    }

    @Override
    int reduce() {
        HashSet<Integer> toDelete = new HashSet<>();
        for(int v : algorithm.getGraph().getVertices()){
            if(algorithm.getGraph().getNeighbors(v).size() < this.threshold)
                toDelete.add(v);
        }
        toDelete.forEach(v -> algorithm.deleteVertex(v));
        return toDelete.size();
    }

    @Override
    boolean useUntilExhaustion() {
        return true;
    }

    @Override
    public String getRuleName() {
        return logName;
    }
}
