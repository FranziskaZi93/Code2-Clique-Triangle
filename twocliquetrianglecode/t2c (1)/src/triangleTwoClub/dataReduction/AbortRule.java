package triangleTwoClub.dataReduction;

import triangleTwoClub.ITriangleTwoClubAlgorithm;

/**
 * An AbortRule identifies NO-instances during the algorithm.
 * If the current graph is a NO-instance then this graph can not find a valid or better solution.
 */
public abstract class AbortRule extends DataReductionRule {

    public AbortRule(ITriangleTwoClubAlgorithm algorithm) {
        super(algorithm);
    }

    /**
     * Checks if the current graph is still valid or otherwise if the branch can get aborted
     * @return true id the graph is still valid
     */
    abstract boolean isValid();

    @Override
    protected void applyRule() {
        if(!isValid()){
            algorithm.abortBranch();
            algorithm.getOutput().log(getRuleName(), 1L);
        }
    }

    @Override
    public boolean applyInitially() {
        return false;
    }
}
