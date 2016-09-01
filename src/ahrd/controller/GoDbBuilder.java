package ahrd.controller;

import java.io.*;
import java.net.URL;
import java.nio.channels.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.HashMap;
//import java.util.Iterator;
import java.util.Map;
//import java.util.Map.Entry;
import java.util.zip.GZIPInputStream;
import org.apache.commons.compress.archivers.tar.*;

import ahrd.model.GOterm;

public class GoDbBuilder {

	public static final String reviewedUniProtURL = "ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz";
	public static final String reviewedUniProtFilePath = "data/uniprot_sprot.dat.gz";
	public static final Pattern reviewedUniProtGoAnnotationRegex = Pattern.compile("^DR\\s{3}GO;\\s+(?<goTerm>GO:\\d{7}).*");
	public static final String geneOntologyMYSQLdumpURL = "http://archive.geneontology.org/latest-termdb/go_daily-termdb-tables.tar.gz";
	public static final String geneOntologyMYSQLdumpFilePath = "data/go_daily-termdb-tables.tar.gz";
	public static final Pattern geneOntologyMYSQLdumpTermTableRegex = Pattern.compile("^(?<id>\\d+)\\t(?<name>[^\\t]+)\\t(?<termtype>[^\\t]+)\\t(?<acc>GO:\\d{7})\\t(?<isobsolete>\\d)\\t(?<isroot>\\d)\\t(?<isrelation>\\d)$");
	public static final Pattern geneOntologyMYSQLdumpGraphPathRegex = Pattern.compile("^(?<id>\\d+)\\t(?<term1id>\\d+)\\t(?<term2id>\\d+)\\t(?<relationshiptypeid>\\d+)\\t(?<distance>\\d+)\\t(?<relationdistance>\\d+)$");
	//public static final Pattern geneOntologyMYSQLdumpGraphPathRegex = Pattern.compile("^(?<id>\\d+)\\t(?<term1id>\\d+)\\t(?<term2id>\\d+).*");
	private static Map<String, GOterm> goDB = new HashMap<String, GOterm>();
	private static Map<String, Integer> annotations = new HashMap<String, Integer>();
	
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
					String line, acc = new String();
					while ((line = termtablebr.readLine()) != null) {
						Matcher m = geneOntologyMYSQLdumpTermTableRegex.matcher(line);
						if (m.matches()) {
							acc = m.group("acc");
							if (goDB.containsKey(acc)) {
								System.out.println("Redundancy in GO term data found.");
							} else {
								GOterm term = new GOterm(acc, m.group("name"), m.group("termtype"));
								goDB.put(acc, term);
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
					while ((line = br.readLine()) != null) {
						linecount++;
						Matcher m = geneOntologyMYSQLdumpGraphPathRegex.matcher(line);
						if (m.matches()) {
							matchcount++;
						}
					}
					System.out.println(matchcount + " of " + linecount + " lines matched");
				}
				entry = tais.getNextTarEntry();
			}
		} finally {
			tais.close();
		}
	}
}
