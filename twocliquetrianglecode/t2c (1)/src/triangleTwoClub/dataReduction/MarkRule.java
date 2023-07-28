package triangleTwoClub.dataReduction;

import triangleTwoClub.ITriangleTwoClubAlgorithm;

/**
 * A MarkRule marks vertices that have to be part of a branch solution based on the already marked vertices.
 */
public abstract class MarkRule extends DataReductionRule {

    public MarkRule(ITriangleTwoClubAlgorithm algorithm) {
        super(algorithm);
    }

    abstract int mark();

    @Override
    protected void applyRule() {
        algorithm.getOutput().log(getRuleName(), mark());
    }

    @Override
    public boolean applyInitially() {
        return false;
    }
}
