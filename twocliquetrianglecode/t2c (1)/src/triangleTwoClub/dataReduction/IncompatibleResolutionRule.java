package triangleTwoClub.dataReduction;

import triangleTwoClub.ITriangleTwoClubAlgorithm;

import java.util.HashSet;

/**
 * Reduction Rule 2 - Incompatible Resolution Rule
 * This DeleteRule removes vertices that are incompatible to marked ones.
 * A marked vertex has to be part of the currents branches solution and therefore we can remove those who are
 * incompatible to marked ones.
 */
public class IncompatibleResolutionRule extends DeleteRule{

    public IncompatibleResolutionRule(ITriangleTwoClubAlgorithm algorithm) {
        super(algorithm);
    }

    @Override
    int reduce() {
        HashSet<Integer> toDelete = new HashSet<>();
        for(int v : algorithm.getMarked()){
            for(int w : algorithm.getGraph().getVertices()){
                if(!algorithm.getConditionChecker().verticesAreCompatible(v, w))
                    toDelete.add(w);
            }
        }
        toDelete.forEach(v -> algorithm.deleteVertex(v));
        return toDelete.size();
    }

    @Override
    public String getRuleName() {
        return "IncompatibleResolutionRule";
    }

    @Override
    boolean useUntilExhaustion() {
        return true;
    }
}
