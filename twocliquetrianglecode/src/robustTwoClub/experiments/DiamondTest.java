package robustTwoClub.experiments;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;


public class DiamondTest {

	public static void main (String[] args) {
		int NODES = 16;// CHANGE THIS TO ALTER NUMBER OF VERTICES IN GENERATED GRAPHS
		int ITER = 5000;	// CHANGE THIS TO ALTER NUMBER OF GRAPHS GENERATED
		boolean OUTPUT = false;	 // Set this to 'true' to generate output files in .dim format
		boolean TYPEDIM = false; // Set to false to let generated files end on .edges instead of .dim
		int maxEdges = 0; int minEdges = NODES * NODES;
		int maxNEdges = 0; int minNEdges = NODES * NODES;
		int maxNNEdges = 0; int minNNEdges = NODES * NODES;
		int maxDiamond = 0; int minDiamond = NODES * NODES;
		int maxMinDegree = 0; int minMinDegree = NODES * NODES;
		int maxMaxDegree = 0; int minMaxDegree = NODES * NODES;
		int maxTriangles = 0; int minTriangles = NODES * NODES;
		for (int i = 0; i < ITER; i++) {
			boolean[][] graph = generate22Club(NODES);
			if (!check22Club(graph)) {
				System.out.println("ERROR: Generated graph is not a 2,2-Club!");
				System.exit(0);
			}
			int edges = countEdges(graph);
			int diamond = largestDiamond(graph);
			int degree[] = minMaxDegree(graph);
			int necEdges = countNecessaryEdges(graph);
			int NNEdges = edges - necEdges;
			int triangles = countTriangles(graph);
			System.out.println("#"+i+" Edges: " + edges + "  Largest diamond: " + diamond + 
					"  MinDeg: " + degree[0] + "  MaxDeg: " + degree[1] + "  Necessary Edges: " + necEdges + 
					"  Non-nec. edges: " + NNEdges + "  Triangles: " + triangles);
			if (edges > maxEdges) maxEdges = edges;
			if (edges < minEdges) minEdges = edges;
			if (necEdges > maxNEdges) maxNEdges = necEdges;
			if (necEdges < minNEdges) minNEdges = necEdges;
			if (NNEdges > maxNNEdges) maxNNEdges = NNEdges;
			if (NNEdges < minNNEdges) minNNEdges = NNEdges;
			if (diamond > maxDiamond) maxDiamond = diamond;
			if (diamond < minDiamond) minDiamond = diamond;
			if (degree[0] > maxMinDegree) maxMinDegree = degree[0];
			if (degree[0] < minMinDegree) minMinDegree = degree[0];
			if (degree[1] > maxMaxDegree) maxMaxDegree = degree[1];
			if (degree[1] < minMaxDegree) minMaxDegree = degree[1];
			if (triangles > maxTriangles) maxTriangles = triangles;
			if (triangles < minTriangles) minTriangles = triangles;
			if (OUTPUT) {
				if (TYPEDIM) createDIM(graph, "graph"+i+".dim");
				else createDIM(graph, "graph"+i+".edges");
			}	
		}
		System.out.println("\nEdge range: " + minEdges + " - " + maxEdges);
		System.out.println("Range of necessary edges: " + minNEdges + " - " + maxNEdges);
		System.out.println("Range of non-necessary edges: " + minNNEdges + " - " + maxNNEdges);
		System.out.println("Largest diamond range: " + minDiamond + " - " + maxDiamond);
		System.out.println("Min degree range: " + minMinDegree + " - " + maxMinDegree);
		System.out.println("Max degree range: " + minMaxDegree + " - " + maxMaxDegree);
		System.out.println("Triangle range: " + minTriangles + " - " + maxTriangles);
	}
	
	private static int countTriangles (boolean[][] graph) {
		int triangles = 0;
		for (int i = 0; i < graph.length; i++)
			for (int j = i+1; j < graph.length; j++) 
				if (graph[i][j] == true)
					for (int k = j+1; k < graph.length; k++) 
						if (graph[i][k] == true)
							if (graph[j][k] == true)
								triangles++;
		return triangles;
	}
	
	private static int countNecessaryEdges (boolean[][] graph) {
		// NOTE: An edge, which is not necessary can be removed without destroying the 2,2-club property,
		// but that doesn't mean that we can remove ALL non-necessary edges at once, just each alone.
		int necEdges = 0;
		for (int i = 0; i < graph.length; i++)
			for (int j = i+1; j < graph.length; j++) 
				if (graph[i][j] == true) {
					boolean edgeNecessary = false;
					// Check whether i and j need this edge
					int neighbors = 0;
					for (int k = 0; k < graph.length; k++)
						if (graph[i][k] == true && graph[j][k] == true)
							neighbors++;
					if (neighbors < 2) {
						necEdges++;
						continue;
					}
					for (int k = 0; k < graph.length; k++)
						if (graph[i][k] == true && k != j) {
							if (graph[j][k] == true) continue;
							neighbors = 0;
							for (int l = 0; l < graph.length; l++)
								if (graph[j][l] == true && graph[k][l] == true && l != i)
									neighbors++;
							if (neighbors < 2) {
								necEdges++;
								edgeNecessary = true;
								break;
							}
						}
					if (!edgeNecessary)
						for (int k = 0; k < graph.length; k++)
							if (graph[j][k] == true && k != i) {
								if (graph[i][k] == true) continue;
								neighbors = 0;
								for (int l = 0; l < graph.length; l++)
									if (graph[i][l] == true && graph[k][l] == true && l != j)
										neighbors++;
								if (neighbors < 2) {
									necEdges++;
									break;
								}
							}
				}
		return necEdges;
	}
	
	private static int[] minMaxDegree (boolean[][] graph) {
		int[] degree = new int[2];
		degree[0] = graph.length;
		degree[1] = 0;
		for (int node = 0; node < graph.length; node++) {
			int neighbors = 0;
			for (int k = 0; k < graph.length; k++)
				if (graph[node][k] == true)
					neighbors++;
			if (neighbors > degree[1])
				degree[1] = neighbors;
			if (neighbors < degree[0])
				degree[0] = neighbors;
		}
		return degree;
	}
	
	private static boolean check22Club (boolean[][] graph) {
		for (int i = 0; i < graph.length; i++)
			for (int j = 0; j < graph.length; j++)
				if (i != j && graph[i][j] == false) {
					int neighbors = 0;
					for (int k = 0; k < graph.length; k++)
						if (graph[i][k] == true && graph[j][k] == true)
							neighbors++;
					if (neighbors < 2) return false;
				}
		return true;
	}
	
	private static int countEdges (boolean[][] graph) {
		int edges = 0;
		for (int i = 0; i < graph.length; i++)
			for (int j = 0; j < graph.length; j++)
				if (graph[i][j] == true) edges++;
		edges /= 2;
		return edges;
	}
	
	private static int largestDiamond (boolean[][] graph) {
		int diamond = 0;
		for (int i = 0; i < graph.length; i++)
			for (int j = i+1; j < graph.length; j++) {
				int neighbors = 0;
				for (int k = 0; k < graph.length; k++)
					if (graph[i][k] == true && graph[j][k] == true)
						neighbors++;
				if (neighbors > 1 && neighbors + 2 > diamond)
					diamond = neighbors + 2;
			}
		return diamond;
	}
	
	private static boolean[][] generate22Club (int n) {
		boolean[][] adjM = new boolean[n][n];	// adjacency matrix
		boolean[][] agrM = new boolean[n][n];	// agreeability matrix
		int nagr = n * (n - 1);					// number of non-agreeable vertex pairs (each counted twice!)
		int[] nagrV = new int[n];				// number of non-agreeable vertices for each vertex
		for (int i = 0; i < n; i++)
			nagrV[i] = n-1;
		while (nagr > 0) {							// while there is some non-agreeable pair
			int pair = (int)(Math.random()*nagr);	// pick one pair uniformly at random
			// Determine picked pair from pair number
			int row = 0; int col = 0;				
			while (nagrV[row] <= pair) {
				pair -= nagrV[row];
				row++;
			}	
			while (pair > 0) {
				if (row == col) {
					col++;		
					continue;	// skip diagonal
				}
				if (agrM[row][col] == false) pair--;
				col++;
			}
			while (agrM[row][col] == true)
				col++;
			if (row == col) col++;
			while (agrM[row][col] == true)
				col++;
			// Insert edge for picked pair
			adjM[row][col] = true;
			adjM[col][row] = true;
			// Mark pair as agreeable
			agrM[row][col] = true;
			nagrV[row]--;
			agrM[col][row] = true;
			nagrV[col]--;
			nagr -= 2;
			// Check for other pairs made agreeable
			for (int i = 0; i < n; i++) {
				if (i != row && adjM[col][i] == true) {	// found neighbor i of col
					if (agrM[row][i] == false) {		// if row and i weren't agreeable before
						boolean commonNeighbor = false;	// check if they are now by looking for common neighbor besides col
						for (int j = 0; j < n; j++) {
							if (j != col && adjM[row][j] == true && adjM[i][j] == true)
								commonNeighbor = true;
						}
						if (commonNeighbor) {			// update agreeability information
							agrM[row][i] = true;
							nagrV[row]--;
							agrM[i][row] = true;
							nagrV[i]--;
							nagr -= 2;
						}
					}
				}	
			}
			// Same procedure for neighbors of row
			for (int i = 0; i < n; i++) {
				if (i != col && adjM[row][i] == true) {	// found neighbor i of row
					if (agrM[col][i] == false) {		// if col and i weren't agreeable before
						boolean commonNeighbor = false;	// check if they are now by looking for common neighbor besides row
						for (int j = 0; j < n; j++) {
							if (j != row && adjM[col][j] == true && adjM[i][j] == true)
								commonNeighbor = true;
						}
						if (commonNeighbor) {			// update agreeability information
							agrM[col][i] = true;
							nagrV[col]--;
							agrM[i][col] = true;
							nagrV[i]--;
							nagr -= 2;
						}
					}
				}	
			}
		}
		return adjM;
	}
	
	private static void createDIM (boolean[][] graph, String filename) {
		File file = new File(filename);
		try {
			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
			for (int i = 0; i < graph.length; i++)
				for (int j = i+1; j < graph.length; j++)
					if (graph[i][j] == true)	
						bw.write(i + " " + j + "\n");
			bw.flush();
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}
