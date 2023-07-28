package robustTwoClub.latex;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Locale;

public class TableMaker {

	private static DecimalFormat formatTime; // = new DecimalFormat("#,##0.00");
	private static DecimalFormat formatDensity; // = new DecimalFormat("0.0000");
	
	private static ArrayList<String> lines;
	private static ArrayList<String> csvLines;
	private static String curLine;
	private static String graphName;
	
	public static void init(String name) {
		final Locale locale = new Locale("en", "US");

		formatTime = (DecimalFormat) NumberFormat.getNumberInstance(locale);
		formatTime.applyPattern("###,###,###.##");

		formatDensity = (DecimalFormat) NumberFormat.getNumberInstance(locale);
		formatDensity.applyPattern("###,###,###.####");
		
		lines = new ArrayList<String>();
		csvLines = new ArrayList<String>();
		lines.add("\\documentclass{article}");
		lines.add("\\begin{document}");
		lines.add("\\begin{figure}");
		lines.add("\\begin{tabular}[c]{|c|ccccc|ccccccc|}");
		lines.add("\\hline");
		lines.add("t & n & m & density & $\\Delta_{s}$ & d & time & time$_{ST}$ & K & K$_{hard}$"); 
		lines.add("& n$_{hard}$ & branches & ops \\\\");
		lines.add("\\hline");
		csvLines.add("t & n & m & density & $\\Delta_{s}$ & d & time & time$_{ST}$ & K & K$_{hard}$ & n$_{hard}$ & branches & ops ");
		curLine = "";
		graphName = name;
	}
	
	public static void addSolutionInfo(int t, int n, int m, double density, int sDelta, int degeneracy) {
		curLine += t + " & ";
		curLine += n + " & ";
		curLine += m + " & ";
		curLine += formatDensity.format(density) + " & ";
		if (sDelta == - 1) curLine += "--- & ";
		else curLine += sDelta + " & ";
		if (degeneracy >= t) curLine += "\\textbf{" + degeneracy + "} & ";
		else curLine += degeneracy + " & ";
	}
	
	public static void addAlgorithmInfo(double time, double btime, int kernels, int hkernels, int hksize,
			int branches, int ops) {
		curLine += formatTime.format(time) + " & ";
		curLine += formatTime.format(btime) + " & ";
		curLine += kernels + " & ";
		curLine += hkernels + " & ";
		curLine += hksize + " & ";
		curLine += branches + " & ";
		curLine += ops;
		csvLines.add(curLine);
		curLine += " \\\\";
		lines.add(curLine);
		curLine = "";
	}
	
	public static void addCliqueLine (int t, int n) {
		lines.add(t+" & " + n + " & " + n*(n-1)/2 + " & 1.00 & & & & & & & & & \\\\");
		csvLines.add(t+" & " + n + " & " + n*(n-1)/2 + " & 1.00 & & & & & & & & &");
	}
	
	public static void finalizeAndOutput(String directory, String filename) {
		lines.add("\\hline");
		lines.add("\\end{tabular}");
		lines.add("\\caption{Results for the graph '"+graphName+"'. "
				+ "If we omit some value of $t$ in a figure, then the solution for this value of $t$ is the same as the solution for the next listed value of $t$. "
				+ "The columns denote the following values: $n$ and $m$ the number of vertices and edges in the solution, "
				+ "'density' the density of the solution, $d$ the largest diamond within the solution, "
				+ "$\\Delta_s$ the number of vertices this solution that were not part of the solution for the previous value of $t$, "
				+ "'time' the time in seconds for sovling the problem for a given value of $t$, "
				+ "'time$_{ST}$' the time spent in the search tree phase, "
				+ "$K$ the number of Turing kernels, $K_{hard}$ the number of hard Turing kernels (where data reduction alone does not suffice), "
				+ "$n_{hard}$ the size of the largest hard kernel, 'branches' the number of search tree nodes and 'ops' is some number of operations (needs more precise description!).}");
		lines.add("\\end{figure}");
		lines.add("\\end{document}");
		File file, file2; File dir;
		if (directory != null) {
			file = new File(directory+"/"+filename+".tex");
			file2 = new File(directory+"/"+filename+".csv");
			dir = new File(directory+"/");
			if (!dir.exists()) dir.mkdir();
		} else {
			file = new File(filename+".tex");
			file2 = new File(filename+".csv");
		}
		try {
			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);
			// Write information to file
			for (String line : lines)
				bw.write(line + "\n");
			bw.flush();
			bw.close();
			FileWriter fw2 = new FileWriter(file2.getAbsoluteFile());
			BufferedWriter bw2 = new BufferedWriter(fw2);
			for (String line : csvLines)
			{
				line = line.replace("&", "\t");
				line = line.replace("$", "");
				line = line.replace("_", "-");
				line = line.replace(" ", "");
				line = line.replace("\\", "");
				line = line.replace("textbf", "");
				line = line.replace("{", "");
				line = line.replace("}", "");
				bw2.write(line + "\n");
			}
			bw2.flush();
			bw2.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(0);
		}
	}
	
}
