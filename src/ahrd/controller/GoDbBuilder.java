package ahrd.controller;

import java.io.*;
import java.net.URL;
import java.nio.channels.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import org.apache.commons.compress.archivers.tar.*;
import java.lang.Math;

import ahrd.model.GOterm;

public class GoDbBuilder {

	public static final String reviewedUniProtURL = "ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz";
	public static final String reviewedUniProtFilePath = "data/uniprot_sprot.dat.gz";
	public static final Pattern reviewedUniProtGoAnnotationRegex = Pattern.compile("^DR\\s{3}GO;\\s+(?<goTerm>GO:\\d{7}).*");
	public static final String geneOntologyMYSQLdumpURL = "http://archive.geneontology.org/latest-termdb/go_daily-termdb-tables.tar.gz";
	public static final String geneOntologyMYSQLdumpFilePath = "data/go_daily-termdb-tables.tar.gz";
	public static final Pattern geneOntologyMYSQLdumpTermTableRegex = Pattern.compile("^(?<id>\\d+)\\t(?<name>[^\\t]+)\\t(?<termtype>[^\\t]+)\\t(?<acc>[^\\t]+)\\t(?<isobsolete>\\d)\\t(?<isroot>\\d)\\t(?<isrelation>\\d)$");
	public static final Pattern geneOntologyMYSQLdumpGraphPathRegex = Pattern.compile("^(?<id>\\d+)\\t(?<term1id>\\d+)\\t(?<term2id>\\d+)\\t(?<relationshiptypeid>1|25)\\t(?<distance>\\d+)\\t(?<relationdistance>\\d+)$");
	private static Map<Integer, GOterm> goDB = new HashMap<Integer, GOterm>();
	private static Map<String, Integer> annotations = new HashMap<String, Integer>();
	private static int ccCount = 0;
	private static int mfCount = 0;
	private static int bpCount = 0;
	
	public static void main(String[] args) throws IOException {
		// Download swiss prot if not already on drive
		if (!new File(reviewedUniProtFilePath).exists()) {
			System.out.println("Downloading " + reviewedUniProtURL); 
			download(reviewedUniProtURL, reviewedUniProtFilePath);
		}
		// Download gene ontology mysql data base dump if not alrady on drive
		if (!new File(geneOntologyMYSQLdumpFilePath).exists()) {
			System.out.println("Downloading " + geneOntologyMYSQLdumpURL);
			download(geneOntologyMYSQLdumpURL, geneOntologyMYSQLdumpFilePath);
		}
		
		//Count swiss prot annotations per GO term
		countSwissProtAnnoations();
		
		// Read term table and fill goDB
		readTermTable();
		
		// Read graph path table and fill the terms ancestry in goDB
		readGraphPath();
		
		// Sum annotation counts to term frequencies
		sumCounts();
		
		// Calculate term probabilities from the term frequencies and then the information content from the probabilities  
		calculateInfoContent();
		
		// for testing
		BufferedWriter out = null;
		System.out.println("Size of GO annotation map: " + annotations.size());
		try {
			out = new BufferedWriter(new FileWriter("data/annotations.csv"));
		    Iterator<Entry<String, Integer>> it = annotations.entrySet().iterator();
		    while (it.hasNext()) {
		        Map.Entry<String, Integer> pair = it.next();
		        out.write(pair.getKey() + "," + pair.getValue() + "\n");

		    }
		    System.out.println("Done!");
		} finally {
			out.close();
		}

	}
	
	private static void download(String url, String file) throws IOException {
		URL website = new URL(url);
		ReadableByteChannel rbc = Channels.newChannel(website.openStream());
		FileOutputStream fos = new FileOutputStream(file);
		fos.getChannel().transferFrom(rbc, 0, Long.MAX_VALUE);
		fos.close();
	}
	
	//Count swiss prot annotations per GO term
	private static void countSwissProtAnnoations() throws FileNotFoundException, IOException {
		String spline, go = new String();
		BufferedReader swissProtIn = null;
		try {
			swissProtIn = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(reviewedUniProtFilePath))));
			while ((spline = swissProtIn.readLine()) != null) {
				Matcher m = reviewedUniProtGoAnnotationRegex.matcher(spline);
				if (m.matches()) {
					go = m.group("goTerm");
					if (annotations.containsKey(go)) {
						annotations.put(go, annotations.get(go) + 1);
					} else {
						annotations.put(go, 1);
					}
				}
			}
		} finally {
			swissProtIn.close();
		}
	}
	
	// Read term table and fill goDB 
	private static void readTermTable() throws FileNotFoundException, IOException {
		TarArchiveInputStream godbtis = null;
		try {
			godbtis = new TarArchiveInputStream(new GZIPInputStream(new FileInputStream(geneOntologyMYSQLdumpFilePath)));
			TarArchiveEntry entry = godbtis.getNextTarEntry();
			while (entry != null) {
				if (entry.getName().equals("go_daily-termdb-tables/term.txt")) {
					BufferedReader termtablebr = null;
					termtablebr = new BufferedReader(new InputStreamReader(godbtis));
					Integer id = 0;
					String line = new String();
					while ((line = termtablebr.readLine()) != null) {
						Matcher m = geneOntologyMYSQLdumpTermTableRegex.matcher(line);
						if (m.matches()) {
							id = Integer.parseInt(m.group("id"));
							if (goDB.containsKey(id)) {
								System.out.println("Redundancy in GO term data found.");
							} else {
								GOterm term = new GOterm(m.group("acc"), m.group("name"), m.group("termtype"));
								goDB.put(id, term);
							}
						}
					}
					System.out.println("Size of goDB:" + goDB.size());
				}
				entry = godbtis.getNextTarEntry();
			}
		} finally {
			godbtis.close();
		}
	}

	// Read graph path table and fill the terms ancestry in goDB
	private static void readGraphPath() throws FileNotFoundException, IOException {
		TarArchiveInputStream tais = null;
		try {
			tais = new TarArchiveInputStream(new GZIPInputStream(new FileInputStream(geneOntologyMYSQLdumpFilePath)));
			TarArchiveEntry entry = tais.getNextTarEntry();
			BufferedReader br = null;			
			while (entry != null) {
				if (entry.getName().equals("go_daily-termdb-tables/graph_path.txt")) {
					br = new BufferedReader(new InputStreamReader(tais));
					String line = new String();
					Integer linecount = 0;
					Integer matchcount = 0;
					Integer failcount = 0;
					while ((line = br.readLine()) != null) {
						linecount++;
						Matcher m = geneOntologyMYSQLdumpGraphPathRegex.matcher(line);
						if (m.matches()) {
							matchcount++;
							GOterm term = goDB.get(Integer.parseInt(m.group("term2id")));
							GOterm parent = goDB.get(Integer.parseInt(m.group("term1id")));
							if (term != null && parent != null) {
								term.addTermToAncestry(parent);	
							} else {
								failcount++;
							}
						}
					}
					System.out.println(matchcount + " of " + linecount + " lines matched");
					System.out.println("Failcount: "+ failcount);
				}
				entry = tais.getNextTarEntry();
			}
		} finally {
			tais.close();
		}
	}
	
	// Sum annotation counts to term frequencies
	private static void sumCounts() {
		Iterator<Map.Entry<Integer, GOterm>> mapIter = goDB.entrySet().iterator();
		GOterm term = null;
		Set<GOterm> ancestry = null;
		String acc = new String();
		int count = 0;
		while (mapIter.hasNext()) {
			Map.Entry<Integer, GOterm> entry = mapIter.next();
			term = entry.getValue();
			acc = term.getAccession();
			if (annotations.containsKey(acc)) {
				count = annotations.get(acc);
				switch(term.getOntology()) {
				case "cellular_component":	ccCount = ccCount + count;
											break;
				case "biological_process":	bpCount = bpCount + count;
											break;
				case "molecular_function":	mfCount = mfCount + count;
											break;
				default: //meta term: Not counted towards any of the three ontologies
				}
				ancestry = term.getAncestry();
				// Because each term is part of its own ancestry, its annotation count does not have to be increase separately
				for(GOterm parent : ancestry) {
					parent.setFrequency(parent.getFrequency() + count);
				}
			}
		}
		System.out.println("ccCount: " + ccCount);
		System.out.println("bpCount: " + bpCount);
		System.out.println("mfCount: " + mfCount);
	}
	
	// Calculate term probabilities from the term frequencies and then the information content from the probabilities
	private static void calculateInfoContent() {
		Iterator<Map.Entry<Integer, GOterm>> mapIter = goDB.entrySet().iterator();
		GOterm term = null;
		while (mapIter.hasNext()) {
			Map.Entry<Integer, GOterm> entry = mapIter.next();
			term = entry.getValue();
			switch(term.getOntology()) {
			case "cellular_component":	term.setProbability((double)term.getFrequency()/ccCount);
										break;
			case "biological_process":	term.setProbability((double)term.getFrequency()/bpCount);
										break;
			case "molecular_function":	term.setProbability((double)term.getFrequency()/mfCount);
										break;
			default: //meta term: Not counted towards any of the three ontologies
			}
			term.setInformationContent(Math.log(term.getProbability()));
		}
	}
}
