package triangleTwoClub.dataReduction;

import triangleTwoClub.ITriangleTwoClubAlgorithm;

/**
 * Reduction Rule 1 - Marked Incompatible Rule
 * This AbortRule checks if all marked vertices in the graph are still compatible to each other.
 * If two marked vertices are not compatible anymore (e.g. as a result of deletions), then this branch can
 * not yield a solution.
 *
 */
public class MarkedIncompatibleRule extends AbortRule{

    public MarkedIncompatibleRule(ITriangleTwoClubAlgorithm algorithm) {
        super(algorithm);
    }

    @Override
    boolean isValid() {
        for(int v : algorithm.getMarked()){
            for(int w : algorithm.getMarked()){
                if(v == w) continue;
                if(!algorithm.getConditionChecker().verticesAreCompatible(v, w))
                    return false;
            }
        }
        return true;
    }

    @Override
    public String getRuleName() {
        return "MarkedIncompatibleRule";
    }
}
