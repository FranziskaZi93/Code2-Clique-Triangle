package robustTwoClub.graph;

public class Edge{
    private Integer v, w;
    public Edge(Integer v, Integer w){
        this.v = v;
        this.w = w;
    }
    public Integer getV(){return  v;}
    public Integer getW(){return  w;}
    public boolean hasNull(){
        return v == null || w == null;
    }
    public Integer getNonNullVertex(){
        return v == null ? w : v;
    }

    @Override
    public boolean equals(Object obj) {
        if(obj == null) return false;
        if(super.equals(obj)) return true;
        if(obj.getClass() == Edge.class) {
            Edge e = (Edge) obj;
            return (e.getV().equals(v) && e.getW().equals(w))
                    || (e.getV().equals(w) && e.getW().equals(v));
        }
        return super.equals(obj);
    }

    public boolean contains(Integer u){
        return v.equals(u) || w.equals(u);
    }

    @Override
    public String toString() {
        return "Edge{" + v + ", " + w + "}";
    }
}