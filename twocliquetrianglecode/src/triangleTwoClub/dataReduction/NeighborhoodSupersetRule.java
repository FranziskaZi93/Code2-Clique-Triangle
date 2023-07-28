package triangleTwoClub.dataReduction;

import robustTwoClub.graph.RtcGraph;
import triangleTwoClub.ITriangleTwoClubAlgorithm;

import java.util.HashSet;

/**
 * Let v, u be vertices.
 * If N[u] \supseteq N_2[v], then delete v.
 */
public class NeighborhoodSupersetRule extends DeleteRule{

    public NeighborhoodSupersetRule(ITriangleTwoClubAlgorithm algorithm) {
        super(algorithm);
    }

    @Override
    int reduce() {
        RtcGraph g = algorithm.getGraph();
        HashSet<Integer> toDelete = new HashSet<>();
        vLoop:
        for(Integer v : g.getVertices()){
            HashSet<Integer> twoNeighborhood = g.getTwoNeighbors(v);
            for(Integer u : twoNeighborhood){
                if(g.getNeighbors(u).containsAll(twoNeighborhood)){
                    toDelete.add(v);
                    continue vLoop;
                }
            }
        }
        toDelete.forEach(v -> algorithm.deleteVertex(v));
        return toDelete.size();
    }

    @Override
    boolean useUntilExhaustion() {
        return true;
    }

    @Override
    public String getRuleName() {
        return "NeighborhoodSupersetRule";
    }
}
