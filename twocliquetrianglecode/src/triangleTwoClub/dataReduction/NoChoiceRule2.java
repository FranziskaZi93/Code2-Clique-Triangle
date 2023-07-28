package triangleTwoClub.dataReduction;

import triangleTwoClub.ITriangleTwoClubAlgorithm;

import java.util.LinkedList;

/**
 * Reduction Rule 6 - No Choice Rule 2
 * If two marked vertices v and w are non-adjacent und share exactly one neighbor u, then mark u.
 */
public class NoChoiceRule2 extends MarkRule{

    public NoChoiceRule2(ITriangleTwoClubAlgorithm algorithm) {
        super(algorithm);
    }

    @Override
    int mark() {
        LinkedList<Integer> toMark = new LinkedList<>();

        Integer u;
        for(int v : algorithm.getMarked()){
            for(int w : algorithm.getMarked()){
                if(!algorithm.isEdgeConditionEnabled() || !algorithm.getGraph().adjacent(v, w)){
                    u = getSingleCommonNeighbor(v, w);
                    if(u != null)
                        toMark.add(u);
                }
            }
        }

        toMark.forEach( v -> algorithm.markVertex(v) );
        return toMark.size();
    }

    @Override
    public String getRuleName() {
        return "NoChoiceRule2";
    }

    private Integer getSingleCommonNeighbor(int v, int w){
        // determine if there is a single common neighbor
        if(algorithm.getConditionChecker().countCommonNeighbors(v, w, 2) != 1)
            return null;

        // find single common neighbor
        for(int u : algorithm.getGraph().getNeighbors(v))
            if(algorithm.getGraph().getNeighbors(w).contains(u))
                return u;

        return null;
    }
}
