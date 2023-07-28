package robustTwoClub.graph;

import java.util.*;

public class Triangle {

    private ArrayList<Integer> vertices;

    public Triangle(int v, int w, int u){
        vertices = new ArrayList<>(3);
        vertices.add(v);
        vertices.add(w);
        vertices.add(u);
        Collections.sort(vertices);
        if(!isValid()){
            throw new RuntimeException("A triangle can not contain the same vertex twice and must have positive IDs!");
        }
    }

    /**
     * Checks if the triangle exists in the given graph. A triangle exists if the graph contains each of
     * the three triangle vertices as well as edges between them.
     * @param graph a graph
     * @return true if the triangle exists
     */
    public boolean exists(RtcGraph graph){
        if(!graph.contains(vertices.get(0)) || !graph.contains(vertices.get(1)) || !graph.contains(vertices.get(2)))
            return false;
        return graph.adjacent(vertices.get(0), vertices.get(1))
                && graph.adjacent(vertices.get(0), vertices.get(2))
                && graph.adjacent(vertices.get(1), vertices.get(2));
    }

    /**
     * Checks if a vertex is part of the triangle
     * @param v the vertex id
     * @return true if v is part of the triangle
     */
    public boolean contains(int v){
        return vertices.contains(v);
    }

    /**
     * Get the list of vertices of a triangle in ascending order.
     * @return an ArrayList<Integer> containing vertex IDs
     */
    public List<Integer> getVertices(){
        return vertices;
    }

    /**
     * returns a unique triangle identifier for a graph of size n
     * @param n the number of vertices of a graph
     * @return a hash value
     */
    public int hash(int n){
        return (vertices.get(0) * n * n) + (vertices.get(1) * n) + vertices.get(2);
    }

    @Override
    public int hashCode(){
        // TODO: use own hash
        return super.hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        try{
            Triangle t = (Triangle) obj;
            List<Integer> v = t.getVertices();
            return vertices.get(0).equals(v.get(0))
                    && vertices.get(1).equals(v.get(1))
                    && vertices.get(2).equals(v.get(2));
        }catch(ClassCastException e){
            return false;
        }
    }

    private boolean isValid(){
        boolean distinctVertices = !vertices.get(0).equals(vertices.get(1))
                && !vertices.get(0).equals(vertices.get(2))
                && !vertices.get(1).equals(vertices.get(2));
        boolean positive = vertices.get(0) >= 0 && vertices.get(1) >= 0 && vertices.get(2) >= 0;
        return distinctVertices && positive;
    }

    /**
     * Computes a list of every distinct triangle for a given graph
     * @param graph the graph
     * @param timeLimit the time limit for computation in ms
     * @return a list of triangles, or null if timeLimit was reached
     */
    public static List<Triangle> getTriangleList(RtcGraph graph, long timeLimit){
        List<Triangle> triangles = new LinkedList<>();
        long startTime = System.currentTimeMillis();
        HashSet<Integer> commonNeighbors;
        Triangle t;
        for(Integer v : graph.getVertices()){
            for(Integer w : graph.getNeighbors(v)){
                if(w <= v) continue;
                if(timeLimit > 0 && System.currentTimeMillis() > startTime + timeLimit)
                    return null;
                commonNeighbors = graph.getCommonNeighbors(v, w);
                for(Integer u : commonNeighbors){
                    if(u <= w) continue;
                    t = new Triangle(v, w, u);
                    triangles.add(t);
                }
            }
        }

        return triangles;
    }

    /**
     * Computes a list of every distinct triangle for a given graph
     * @param graph the graph
     * @return a list of triangles, or null if timeLimit was reached
     */
    public static List<Triangle> getTriangleList(RtcGraph graph){
        return getTriangleList(graph, -1);
    }

    /**
     * Computes the triangles of a given graph. A triangle of vertices v, w, u will be present in the HashMap under
     * each key v, w, u. The vertex IDs are used as keys.
     * @param graph the graph
     * @return a HashMap with triangles of every vertex
     */
    public static HashMap<Integer, HashSet<Triangle>> getTriangleHashMap(RtcGraph graph){
        return getTriangleHashMap(graph, -1);
    }

    /**
     * Computes the triangles of a given graph. A triangle of vertices v, w, u will be present in the HashMap under
     * each key v, w, u. The vertex IDs are used as keys.
     * @param graph the graph
     * @param timeLimit the time limit for computation in ms
     * @return a HashMap with triangles of every vertex, or null if the timeLimit was reached
     */
    public static HashMap<Integer, HashSet<Triangle>> getTriangleHashMap(RtcGraph graph, long timeLimit){
        List<Triangle> list = Triangle.getTriangleList(graph, timeLimit);
        if(list == null)
            return null;
        HashMap<Integer, HashSet<Triangle>> map = new HashMap<>();

        for(Triangle t : list){
            for(Integer v : t.getVertices()){
                if(!map.containsKey(v))
                    map.put(v, new HashSet<>());
                map.get(v).add(t);
            }
        }

        return map;
    }

    public static HashMap<Integer, HashSet<Triangle>> getLocalTriangleHashMap(RtcGraph graph,
                                                                              HashSet<Integer> neighborhood){
        HashMap<Integer, HashSet<Triangle>> triangles = new HashMap<>();
        for(int v : neighborhood){
            for(int w : neighborhood){
                if(v <= w) continue;
                for(int u : graph.getCommonNeighbors(v, w)){
                    if(u > w && neighborhood.contains(u)){
                        Triangle t = new Triangle(v, w, u);
                        if(!triangles.containsKey(v))
                            triangles.put(v, new HashSet<>());
                        triangles.get(v).add(t);
                        if(!triangles.containsKey(w))
                            triangles.put(w, new HashSet<>());
                        triangles.get(w).add(t);
                        if(!triangles.containsKey(u))
                            triangles.put(u, new HashSet<>());
                        triangles.get(u).add(t);
                    }
                }
            }
        }
        return triangles;
    }
}

