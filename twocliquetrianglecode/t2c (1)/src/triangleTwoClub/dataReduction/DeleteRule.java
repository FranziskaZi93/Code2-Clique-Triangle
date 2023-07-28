package triangleTwoClub.dataReduction;

import triangleTwoClub.ITriangleTwoClubAlgorithm;

/**
 * A DeleteRule removes vertices from the graph that can not be part of a solution.
 */
public abstract class DeleteRule extends DataReductionRule {

    public DeleteRule(ITriangleTwoClubAlgorithm algorithm) {
        super(algorithm);
    }

    /**
     * Applies the deletion rule.
     * @return the number of vertices deleted.
     */
    abstract int reduce();

    /**
     * A DeleteRule can be applied once or used until exhaustion.
     * @return true if the rule shall be allied until exhaustion
     */
    abstract boolean useUntilExhaustion();

    /**
     * Uses the reduction rule until exhaustion.
     * @return the total number of deletions.
     */
    private int exhaust(){
        int total = 0, current;
        do{
            current = reduce();
            total += current;
        }while(current > 0);
        return total;
    }

    @Override
    protected void applyRule() {
        if(useUntilExhaustion())
            algorithm.getOutput().log(getRuleName(), exhaust());
        else
            algorithm.getOutput().log(getRuleName(), reduce());
    }

    @Override
    public boolean applyInitially() {
        return true;
    }
}
