package triangleTwoClub.dataReduction;

import robustTwoClub.graph.RtcGraph;
import triangleTwoClub.ITriangleTwoClubAlgorithm;

import java.util.HashSet;

/**
 * Reduction Rule 7 - Maximum Matching Rule
 * Computed the size of a maximum matching for the conflict graph. If the maximum matching is too small, then this
 * branch can not yield a better solution than the current one.
 */
public class MaximumMatchingRule extends AbortRule {

    public MaximumMatchingRule(ITriangleTwoClubAlgorithm algorithm) {
        super(algorithm);
    }

    @Override
    boolean isValid() {
        int limit = algorithm.getGraph().size() - algorithm.getBestSolutionSize() - 1;
        int matchingSize = matchingBoundVertexCover(algorithm.getGraph(), algorithm.getConflictGraph(), limit);
        return !(algorithm.getGraph().size() - matchingSize <= algorithm.getBestSolutionSize());
    }

    @Override
    public String getRuleName() {
        return "MaximumMatchingRule";
    }

    /**
     * Uses the matching bound approximation for vertex cover.
     * @param working the working graph
     * @param conflict the conflict graph
     * @return the minimum number of vertices you need to delete for a triangle 2 club
     */
    public static int matchingBoundVertexCover(RtcGraph working, RtcGraph conflict, int limit){
        // 2-approximation for vertex cover of conflict graph
        // Minimal vertex cover of conflict graph is lower bound on vertices to delete
        HashSet<Integer> marks = new HashSet<Integer>();
        int size = 0;
        for (int v : conflict.getVertices())
            if (working.contains(v) && !marks.contains(v))
                for (int w : conflict.getNeighbors(v))
                    if (working.contains(w) && w > v)					// every remaining conflict is an edge
                        if (!marks.contains(v) && !marks.contains(w)) {	// if the edge is not covered
                            marks.add(v); marks.add(w); size++;			// cover edge with both end points
                            if (size > limit) return size;				// if current size exceeds allowed limit
                            break;										// v us marked, so the condition of the if-statement is never met again.
                        }												// of deletions we do not need to carry on
        return size;
    }
}
