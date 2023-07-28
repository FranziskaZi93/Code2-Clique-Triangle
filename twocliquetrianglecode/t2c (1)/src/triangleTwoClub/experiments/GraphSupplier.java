package triangleTwoClub.experiments;

import robustTwoClub.graph.RtcGraph;

import java.util.ArrayList;
import java.util.List;

public class GraphSupplier {

    private static String dataRepositoryPath =
            "C:\\Users\\phili\\Documents\\Uni\\Unterlagen\\Bachelorarbeit 2-Clubs\\data\\";
    private static String otherGraphs = "C:\\Users\\phili\\Documents\\Uni\\Job Triangl2Club\\graphs\\";


    public static List<Pair<String, RtcGraph>> getPaperGraphs(){
        ArrayList<Pair<String, RtcGraph>> graphs = new ArrayList<>();
        RtcGraph g;

        // Karate
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/karate.graph", 1);
        assertEqualsInt(34, g.size());
        assertEqualsInt(78, g.getEdgeCount());
        graphs.add(new Pair<>("Karate", g));

        // Dolphins
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/dolphins.graph", 1);
        assertEqualsInt(62, g.size());
        assertEqualsInt(159, g.getEdgeCount());
        graphs.add(new Pair<>("Dolphins", g));

        // Lesmis
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/lesmis.graph", 1);
        assertEqualsInt(77, g.size());
        assertEqualsInt(254, g.getEdgeCount());
        graphs.add(new Pair<>("Lesmis", g));

        // Polbooks
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/polbooks.graph", 1);
        assertEqualsInt(105, g.size());
        assertEqualsInt(441, g.getEdgeCount());
        graphs.add(new Pair<>("Polbooks", g));

        // Adjnoun
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/adjnoun.graph", 1);
        assertEqualsInt(112, g.size());
        assertEqualsInt(425, g.getEdgeCount());
        graphs.add(new Pair<>("Adjnoun", g));

        // Football
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/football.graph", 1);
        assertEqualsInt(115, g.size());
        assertEqualsInt(613, g.getEdgeCount());
        graphs.add(new Pair<>("Football", g));

        // Jazz
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/jazz.graph", 1);
        assertEqualsInt(198, g.size());
        assertEqualsInt(2742, g.getEdgeCount());
        graphs.add(new Pair<>("Jazz", g));

        // Huck
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/huck.col", 2);
        assertEqualsInt(74, g.size());
        assertEqualsInt(301, g.getEdgeCount());
        graphs.add(new Pair<>("Huck", g));

        // Jean
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/jean.col", 2);
        assertEqualsInt(80, g.size());
        assertEqualsInt(254, g.getEdgeCount());
        graphs.add(new Pair<>("Jean", g));

        // 3-FullIns_3
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/3-FullIns_3.col", 2);
        assertEqualsInt(80, g.size());
        assertEqualsInt(346, g.getEdgeCount());
        graphs.add(new Pair<>("3-FullIns_3", g));

        // David
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/david.col", 2);
        assertEqualsInt(87, g.size());
        assertEqualsInt(406, g.getEdgeCount());
        graphs.add(new Pair<>("David", g));

        // Mug88_1
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/mug88_1.col", 2);
        assertEqualsInt(88, g.size());
        assertEqualsInt(146, g.getEdgeCount());
        graphs.add(new Pair<>("Mug88_1", g));

        // 1-FullIns_4
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/1-FullIns_4.col", 2);
        assertEqualsInt(93, g.size());
        assertEqualsInt(593, g.getEdgeCount());
        graphs.add(new Pair<>("1-FullIns_4", g));

        // Mug100_1
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/mug100_1.col", 2);
        assertEqualsInt(100, g.size());
        assertEqualsInt(166, g.getEdgeCount());
        graphs.add(new Pair<>("Mug100_1", g));

        // Mug100_25
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/mug100_25.col", 2);
        assertEqualsInt(100, g.size());
        assertEqualsInt(166, g.getEdgeCount());
        graphs.add(new Pair<>("Mug100_25", g));

        // 4-FullIns_3
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/4-FullIns_3.col", 2);
        assertEqualsInt(114, g.size());
        assertEqualsInt(541, g.getEdgeCount());
        graphs.add(new Pair<>("4-FullIns_3", g));

        // Games120
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/games120.col", 2);
        assertEqualsInt(120, g.size());
        assertEqualsInt(638, g.getEdgeCount());
        graphs.add(new Pair<>("Games120", g));

        // Miles500
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/miles500.col", 2);
        assertEqualsInt(128, g.size());
        assertEqualsInt(1170, g.getEdgeCount());
        graphs.add(new Pair<>("Miles500", g));

        // Anna
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/anna.col", 2);
        assertEqualsInt(138, g.size());
        assertEqualsInt(493, g.getEdgeCount());
        graphs.add(new Pair<>("Anna", g));

        // 5-FullIns_3
        g = new RtcGraph("testGraphs/tkclub_paper_graphs/5-FullIns_3.col", 2);
        assertEqualsInt(154, g.size());
        assertEqualsInt(792, g.getEdgeCount());
        graphs.add(new Pair<>("5-FullIns_3", g));

        return graphs;
    }

    public static List<GraphSupplier.Pair<String, Integer>> getBachelorThesisGraphs(){
        ArrayList<GraphSupplier.Pair<String, Integer>> repoGraphs = new ArrayList<>();

        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/1-FullIns_4.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/3-FullIns_3.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/4-FullIns_3.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/5-FullIns_3.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/adjnoun.graph", 1));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/anna.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/david.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/dolphins.graph", 1));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/football.graph", 1));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/games120.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/huck.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/jazz.graph", 1));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/jean.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/karate.graph", 1));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/lesmis.graph", 1));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/miles500.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/mug100_1.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/mug100_25.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/mug88_1.col", 2));
        repoGraphs.add(new Pair<>("testGraphs/tkclub_paper_graphs/polbooks.graph", 1));
        repoGraphs.add(new Pair<>(otherGraphs + "add20.mtx", 0));
        repoGraphs.add(new Pair<>(otherGraphs + "add32.mtx", 0));
        repoGraphs.add(new Pair<>(otherGraphs + "AuthorsFPruned3.txt", 0));
        repoGraphs.add(new Pair<>(otherGraphs + "celegansneural.graph", 1));
        repoGraphs.add(new Pair<>(otherGraphs + "hep-th.graph", 1));
        repoGraphs.add(new Pair<>(otherGraphs + "netscience.graph", 1));
        repoGraphs.add(new Pair<>(otherGraphs + "PGPgiantcompo.graph", 1));
        repoGraphs.add(new Pair<>(otherGraphs + "polblogs.graph", 1));
        repoGraphs.add(new Pair<>(otherGraphs + "power.graph", 1));

        // data/konect
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "konect\\undirected-simple-small\\adjnoun_adjacency\\out.adjnoun_adjacency_adjacency", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "konect\\undirected-simple-small\\arenas-email\\out.arenas-email", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "konect\\undirected-simple-small\\arenas-jazz\\out.arenas-jazz", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "konect\\undirected-simple-small\\contiguous-usa\\out.contiguous-usa", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "konect\\undirected-simple-small\\dolphins\\out.dolphins", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "konect\\undirected-simple-small\\moreno_zebra\\out.moreno_zebra_zebra", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "konect\\undirected-simple-small\\ucidata-zachary\\out.ucidata-zachary", 0));

        // data/dimacs
        repoGraphs.add(new Pair<>(dataRepositoryPath + "dimacs\\citation_networks\\citationCiteseer.graph", 1));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "dimacs\\citation_networks\\coAuthorsCiteseer.graph", 1));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "dimacs\\citation_networks\\coAuthorsDBLP.graph", 1));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "dimacs\\citation_networks\\coPapersCiteseer.graph", 1));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "dimacs\\citation_networks\\coPapersDBLP.graph", 1));

        // data/network repository
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\bio\\bio-celegans-dir.edges", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\bio\\bio-celegans.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\bio\\bio-diseasome.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\bio\\bio-yeast-protein-inter.edges", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\bio\\bio-yeast.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\brain\\bn-cat-mixed-species_brain_1.edges", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\brain\\bn-fly-drosophila_medulla_1.edges", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\brain\\bn-mouse_retina_1.edges", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\collaboration\\ca-AstroPh.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\collaboration\\ca-CondMat.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\collaboration\\ca-CSphd.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\collaboration\\ca-Erdos992.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\collaboration\\ca-GrQc.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\collaboration\\ca-HepPh.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\collaboration\\ca-netscience.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\collaboration\\ca-sandi_auths.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\econ\\econ-beause.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\econ\\econ-orani678.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\infrastructure\\inf-euroroad.edges", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\infrastructure\\inf-openflights.edges", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\infrastructure\\inf-power.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\infrastructure\\inf-USAir97.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\misc\\comsol.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\misc\\heart2.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\misc\\psmigr_1.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "network repository\\misc\\robot24c1_mat5.mtx", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "STC\\AmazonBooks\\BooksEdges.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "STC\\Actors\\graphAllActors.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "STC\\karateClub\\karateClub.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath + "STC\\LesMis\\LesmisEdges.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "network repository\\social\\soc-academia\\soc-academia.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "network repository\\social\\soc-advogato\\soc-advogato.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "network repository\\social\\soc-anybeat\\soc-anybeat.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "network repository\\social\\soc-brightkite\\soc-brightkite.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "network repository\\social\\soc-buzznet\\soc-buzznet.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "network repository\\social\\soc-delicious\\soc-delicious.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "network repository\\social\\soc-firm-hi-tech\\soc-firm-hi-tech.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "network repository\\social\\soc-hamsterster\\soc-hamsterster.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "network repository\\social\\soc-twitter-follows\\soc-twitter-follows.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "network repository\\social\\soc-wiki-Vote\\soc-wiki-Vote.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "network repository\\social\\soc-youtube-snap\\soc-youtube-snap.txt", 0));
        repoGraphs.add(new Pair<>(dataRepositoryPath +
                "network repository\\social\\tech-RL-caida\\tech-RL-caida.txt", 0));
        return repoGraphs;
    }

    public static List<Pair<String, Integer>> getDataRepoGraphPaths(){
        ArrayList<Pair<String, Integer>> repoGraphs = new ArrayList<>();
        // TODO: load graphs
        return repoGraphs;
    }

    public static List<Pair<String, RtcGraph>> getAllGraphs(){
        ArrayList<Pair<String, RtcGraph>> graphs = (ArrayList<Pair<String, RtcGraph>>) getPaperGraphs();
        RtcGraph g;

        // 3line
        g = new RtcGraph("testGraphs/3line.dimacs", 0);
        graphs.add(new Pair<>("3line", g));

        // 7ttclubTest
        g = new RtcGraph("testGraphs/7ttclubTest.dimacs", 0);
        graphs.add(new Pair<>("7ttclubTest", g));

        // C125-9
        g = new RtcGraph("testGraphs/C125-9.mtx", 0);
        graphs.add(new Pair<>("C125-9", g));

        // c-fat200-1
        g = new RtcGraph("testGraphs/c-fat200-1.mtx", 0);
        graphs.add(new Pair<>("c-fat200-1", g));

        // chesapeake
        g = new RtcGraph("testGraphs/chesapeake.mtx", 0);
        graphs.add(new Pair<>("chesapeake", g));

        // citationCiteseer
        g = new RtcGraph("testGraphs/citationCiteseer.graph", 1);
        graphs.add(new Pair<>("citationCiteseer", g));

        // coPapersDBLP
        g = new RtcGraph("testGraphs/coPapersDBLP.graph", 1);
        graphs.add(new Pair<>("coPapersDBLP", g));

        // delaunay_n10
        g = new RtcGraph("testGraphs/delaunay_n10.mtx", 0);
        graphs.add(new Pair<>("delaunay_n10", g));

        // delaunay_n15
        g = new RtcGraph("testGraphs/delaunay_n15.mtx", 0);
        graphs.add(new Pair<>("delaunay_n15", g));

        // inf-luxembourg_osm
        g = new RtcGraph("testGraphs/inf-luxembourg_osm.mtx", 0);
        graphs.add(new Pair<>("inf-luxembourg_osm", g));

        // MANN-a81
        g = new RtcGraph("testGraphs/MANN-a81.mtx", 0);
        graphs.add(new Pair<>("MANN-a81", g));

        return graphs;
    }

    private static void assertEqualsInt(int a, int b){
        if(a != b)
            throw new RuntimeException("Assertion failed.");
    }

    static void checkPaths(){
        for(Pair<String, Integer> p : getBachelorThesisGraphs()){
            RtcGraph g = new RtcGraph(p.getKey(), p.getValue());
            System.out.println(LowerBoundComparison.getName(p.getKey() + " " + g.size() + " "
                    + g.getEdgeCount()));
        }
    }

    public static class Pair<T1, T2>{
        private T1 key;
        private T2 value;

        public Pair(T1 key, T2 value){
            this.key = key;
            this.value = value;
        }

        public T1 getKey(){
            return key;
        }

        public T2 getValue(){
            return value;
        }
    }
}
