package triangleTwoClub.dataReduction;

import robustTwoClub.graph.Triangle;
import triangleTwoClub.ITriangleTwoClubAlgorithm;

import java.util.HashSet;

/**
 * Modified Reduction Rule 5 - No Choice Rule 3
 *
 * Let v be a vertex and x_v the number of triangles v is part of.
 * If u is a neighbor of v and x_u the number of triangles v and u are part of,
 * then mark u if x_v - x_u < l.
 *
 */
public class NoChoiceRule3 extends MarkRule{

    public NoChoiceRule3(ITriangleTwoClubAlgorithm algorithm) {
        super(algorithm);
    }

    @Override
    int mark() {
        HashSet<Integer> toMark = new HashSet<>();

        int x_v, x_u;
        for(int v : algorithm.getMarked()){
            for(int u : algorithm.getGraph().getNeighbors(v)){
                // skip u if u is already marked
                if(algorithm.getMarked().contains(u))
                    continue;

                // count triangles
                x_v = algorithm.getConditionChecker().countTriangles(v);
                x_u = 0;
                if(algorithm.getTriangles().containsKey(u)){
                    for(Triangle t : algorithm.getTriangles().get(u))
                        if(t.contains(v) && t.exists(algorithm.getGraph()))
                            x_u++;
                }
                if(x_v - x_u < algorithm.getL())
                    toMark.add(u);
            }
        }

        toMark.forEach( v -> algorithm.markVertex(v));
        return toMark.size();
    }

    @Override
    public String getRuleName() {
        return "NoChoiceRule3";
    }

}
