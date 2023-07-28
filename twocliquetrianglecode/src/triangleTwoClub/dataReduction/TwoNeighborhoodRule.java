package triangleTwoClub.dataReduction;

import robustTwoClub.algorithms.TriangleTwoClubAlgorithm;
import triangleTwoClub.ITriangleTwoClubAlgorithm;
import triangleTwoClub.Output;

import java.util.HashSet;
import java.lang.Math;

/**
 *  A Data Reduction Rule to make the graph smaller before the actual algorithm starts.
 *  Based on the value for L, vertices with a degree smaller than a certain threshold get deleted.
 *  This rule does not require data structures provided by the algorithm.
 */
public class TwoNeighborhoodRule extends DeleteRule {
    private int lb;

    public TwoNeighborhoodRule(ITriangleTwoClubAlgorithm algorithm, int lb) {
        super(algorithm);
	    this.lb = lb;
    }

    public TwoNeighborhoodRule(ITriangleTwoClubAlgorithm algorithm){
        this(algorithm,0);
    }

    @Override
    int reduce() {
        HashSet<Integer> toDelete = new HashSet<>();
        for(int v : algorithm.getGraph().getVertices()){
	    /*int degBound= 1;
	    for (int w : algorithm.getGraph().getNeighbors(v)){
		degBound += algorithm.getGraph().degree(w);
	    }
	    if(degBound <= lb){ 
                toDelete.add(v);
		continue;
		}*/
	    if(algorithm.getGraph().sizeOfTwoNeighborhood(v,true) <= this.lb) 
                toDelete.add(v);
        }
        toDelete.forEach(v -> algorithm.deleteVertex(v));
	    //System.out.println("Done with one pass of N2Rule");
	    return toDelete.size();
    }

    @Override
    boolean useUntilExhaustion() {
        return true;
    }

    @Override
    public boolean applyInitially() {
        return false;
    }

    @Override
    public String getRuleName() {
        return "TwoNeighborhoodRule";
    }
}
