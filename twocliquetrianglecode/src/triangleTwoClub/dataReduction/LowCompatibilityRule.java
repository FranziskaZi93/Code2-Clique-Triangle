package triangleTwoClub.dataReduction;

import triangleTwoClub.ITriangleTwoClubAlgorithm;
import triangleTwoClub.Output;

import java.util.HashSet;

/**
 * Reduction Rule 4 - Low Compatibility Rule
 * Remove vertices whose number of compatible vertices is smaller or equal than the size of our current best
 * solution size.
 */
public class LowCompatibilityRule extends DeleteRule{

    private String abortCountName;
    private boolean isAbortCountRegistered;

    public LowCompatibilityRule(ITriangleTwoClubAlgorithm algorithm){
        super(algorithm);
        abortCountName = getRuleName()+"Aborts";
        isAbortCountRegistered = false;
    }

    @Override
    int reduce() {
        int k = algorithm.getBestSolutionSize() + 1;
        HashSet<Integer> toDelete = new HashSet<>();

        for(int v : algorithm.getGraph().getVertices()){
            if(algorithm.getConditionChecker().countCompatibles(v, k, true) < k){
                if(algorithm.getMarked().contains(v)){ // if v is marked, then this branch can get aborted
                    if(!isAbortCountRegistered) registerAbortCount();
                    algorithm.abortBranch();
                    algorithm.getOutput().log(abortCountName, 1);
                    break;
                }
                toDelete.add(v);
            }
        }

        toDelete.forEach(v -> algorithm.deleteVertex(v));
        return toDelete.size();
    }

    @Override
    public String getRuleName() {
        return "LowCompatibilityRule";
    }

    @Override
    boolean useUntilExhaustion() {
        return true;
    }

    @Override
    public boolean applyInitially() {
        return false;
    }

    private void registerAbortCount(){
        algorithm.getOutput().log(abortCountName, 0);
        algorithm.getOutput().setFormat(abortCountName, Output.DataFormat.INTEGER);
        isAbortCountRegistered = true;
    }
}
