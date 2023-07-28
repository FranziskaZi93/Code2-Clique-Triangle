package triangleTwoClub.dataReduction;

import triangleTwoClub.ITriangleTwoClubAlgorithm;

import java.util.HashSet;
import java.util.List;

/**
 * DEPRECATED !!!
 * This rule is a special case of NoChoiceRule 3 for l=1. Use the new one instead.
 *
 * Reduction Rule 5 - No Choice Rule 1
 * If a vertex v is marked and every triangle that v is part of also contains another vertex w, then mark w.
 */
@Deprecated
public class NoChoiceRule extends MarkRule{

    public NoChoiceRule(ITriangleTwoClubAlgorithm algorithm) {
        super(algorithm);
    }

    @Override
    int mark() {
        HashSet<Integer> toMark = new HashSet<>();
        List<Integer> sharedTriangleVertices;
        for(int v : algorithm.getMarked()){
            sharedTriangleVertices = algorithm.getConditionChecker().getAllTriangleSharedVertices(v);
            sharedTriangleVertices.forEach( w -> {
                if(!algorithm.getMarked().contains(w))
                    toMark.add(w);
            });
        }

        toMark.forEach( v -> algorithm.markVertex(v));
        return toMark.size();
    }

    @Override
    public String getRuleName() {
        return "NoChoiceRule1";
    }
}
