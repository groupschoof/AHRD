package ahrd.controller;

import java.io.*;
import java.net.URL;
import java.nio.channels.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import org.apache.commons.compress.archivers.tar.*;
import java.lang.Math;

import ahrd.model.GOterm;

public class GoDbBuilder {

	private static final String reviewedUniProtURL = "ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz";
	private static final String reviewedUniProtFilePath = "data/uniprot_sprot.dat.gz";
	private static final Pattern reviewedUniProtGoAnnotationRegex = Pattern.compile("^DR\\s{3}GO;\\s+(?<goTerm>GO:\\d{7}).*");
	private static final String geneOntologyMYSQLdumpURL = "http://archive.geneontology.org/latest-termdb/go_daily-termdb-tables.tar.gz";
	private static final String geneOntologyMYSQLdumpFilePath = "data/go_daily-termdb-tables.tar.gz";
	private static final Pattern geneOntologyMYSQLdumpTermTableRegex = Pattern.compile("^(?<id>\\d+)\\t(?<name>[^\\t]+)\\t(?<termtype>[^\\t]+)\\t(?<acc>[^\\t]+)\\t(?<isobsolete>\\d)\\t(?<isroot>\\d)\\t(?<isrelation>\\d)$");
	private static final Pattern geneOntologyMYSQLdumpGraphPathRegex = Pattern.compile("^(?<id>\\d+)\\t(?<term1id>\\d+)\\t(?<term2id>\\d+)\\t(?<relationshiptypeid>1|25)\\t(?<distance>\\d+)\\t(?<relationdistance>\\d+)$");
	private static final Pattern geneOntologyMYSQLdumpTermSynonymRegex = Pattern.compile("^(?<termid>\\d+)\\t(?<termsynonym>[^\\t]+)\\t(?<accsynonym>GO:\\d{7})\\t(?<synonymtypeid>[^\\t]+)\\t(?<synonymcategoryid>[^\\t]+)$");
	private static final String bpRootAcc = "GO:0008150";
	private static final String ccRootAcc = "GO:0005575";
	private static final String mfRootAcc = "GO:0003674";
	private static Map<String, Integer> annotations = new HashMap<String, Integer>();
	private static Map<Integer, GOterm> idGoDb = new HashMap<Integer, GOterm>();
	private static Map<String, GOterm> accGoDb = new HashMap<String, GOterm>();
	private static int bpCount = 0;
	private static int ccCount = 0;
	private static int mfCount = 0;
	
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
		
		// Count swiss prot annotations per GO term
		countSwissProtAnnoations();
		
		// Read term table and fill goDB
		readTermTable();
		
		// Read graph path table and fill the terms ancestry in goDB
		readGraphPath();
		
		// Sum annotation counts to term frequencies
		sumCounts();
		
		// Calculate term probabilities from the term frequencies
		// Calculate the information content from the term probabilities
		// Store terms in accession based GO database  
		calculateInfoContent();
		
		// Read in term synonyms and add them to accession based GO database
		readTermSynonyms();
		
		System.out.println("Size of accGoDb: " + accGoDb.size());
		
		// Check consistency of root term annotation count and information content
		checkConsistency();
		
		// for testing
//		BufferedWriter out = null;
//		System.out.println("Size of GO annotation map: " + annotations.size());
//		try {
//			out = new BufferedWriter(new FileWriter("data/annotations.csv"));
//		    Iterator<Entry<String, Integer>> it = annotations.entrySet().iterator();
//		    while (it.hasNext()) {
//		        Map.Entry<String, Integer> pair = it.next();
//		        out.write(pair.getKey() + "," + pair.getValue() + "\n");
//
//		    }
//		} finally {
//			out.close();
//		}
		System.out.println("Done!");
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
							if (idGoDb.containsKey(id)) {
								System.out.println("Redundancy in GO term data found.");
							} else {
								if (m.group("isobsolete").equals("1")) {
									idGoDb.put(id, new GOterm(m.group("acc"), m.group("name"), m.group("termtype"), true));
								} else {
									idGoDb.put(id, new GOterm(m.group("acc"), m.group("name"), m.group("termtype")));
								}
							}
						}
					}
					System.out.println("Size of idGoDb: " + idGoDb.size());
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
							GOterm term = idGoDb.get(Integer.parseInt(m.group("term2id")));
							GOterm parent = idGoDb.get(Integer.parseInt(m.group("term1id")));
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
		Iterator<Map.Entry<Integer, GOterm>> mapIter = idGoDb.entrySet().iterator();
		GOterm term = null;
		Set<GOterm> ancestry = null;
		String acc = new String();
		int count = 0;
		int missingRootCount = 0;
		while (mapIter.hasNext()) {
			Map.Entry<Integer, GOterm> entry = mapIter.next();
			term = entry.getValue();
			if (!term.getObsolete()) {
				acc = term.getAccession();
				if (annotations.containsKey(acc)) {
					count = annotations.get(acc);
					switch(term.getOntology()) {
					case "biological_process":	bpCount = bpCount + count;
												break;
					case "cellular_component":	ccCount = ccCount + count;
												break;
					case "molecular_function":	mfCount = mfCount + count;
												break;
					default: 					break; //meta term: Not counted towards any of the three ontologies
					}
					ancestry = term.getAncestry();
					// Because each term is part of its own ancestry, its annotation count does not have to be increase separately
					boolean missesRoot = true;
					for(GOterm parent : ancestry) {
						parent.setFrequency(parent.getFrequency() + count);
						if (parent.getAccession().equals(bpRootAcc) | parent.getAccession().equals(ccRootAcc) | parent.getAccession().equals(mfRootAcc)) {
							missesRoot = false;
						}
					}
					if (missesRoot) {
						System.out.println("Missing root term in ancestry of: " + term.getAccession());
						missingRootCount++;
					}
				}
			}
		}
		System.out.println("bpCount: " + bpCount);
		System.out.println("ccCount: " + ccCount);
		System.out.println("mfCount: " + mfCount);
		System.out.println("missingRootCount: " + missingRootCount);
	}
	
	// Calculate term probabilities from the term frequencies
	// Calculate the information content from the term probabilities
	// Store terms in accession based GO database
	private static void calculateInfoContent() {
		Iterator<Map.Entry<Integer, GOterm>> mapIter = idGoDb.entrySet().iterator();
		GOterm term = null;
		int metacount = 0;
		while (mapIter.hasNext()) {
			Map.Entry<Integer, GOterm> entry = mapIter.next();
			term = entry.getValue();
			switch(term.getOntology()) {
			case "biological_process":	term.setProbability((double)term.getFrequency()/bpCount);
										term.setInformationContent(-1*Math.log(term.getProbability()));
										accGoDb.put(term.getAccession(), term);
										break;
			case "cellular_component":	term.setProbability((double)term.getFrequency()/ccCount);
										term.setInformationContent(-1*Math.log(term.getProbability()));
										accGoDb.put(term.getAccession(), term);
										break;
			case "molecular_function":	term.setProbability((double)term.getFrequency()/mfCount);
										term.setInformationContent(-1*Math.log(term.getProbability()));
										accGoDb.put(term.getAccession(), term);
										break;
			default: 					metacount++;//meta term: Not counted towards any of the three ontologies; omit storing in accession based term database  
										break;
			}
			
		}
		System.out.println("metacount: " + metacount);
		System.out.println("Finished calculating information content");
	}
	
	// Read in term synonyms and add them to accession based GO database
	private static void readTermSynonyms() throws FileNotFoundException, IOException {
		TarArchiveInputStream tais = null;
		try {
			tais = new TarArchiveInputStream(new GZIPInputStream(new FileInputStream(geneOntologyMYSQLdumpFilePath)));
			TarArchiveEntry entry = tais.getNextTarEntry();
			BufferedReader br = null;			
			while (entry != null) {
				if (entry.getName().equals("go_daily-termdb-tables/term_synonym.txt")) {
					br = new BufferedReader(new InputStreamReader(tais));
					String line = new String();
					Integer linecount = 0;
					Integer matchcount = 0;
					while ((line = br.readLine()) != null) {
						linecount++;
						Matcher m = geneOntologyMYSQLdumpTermSynonymRegex.matcher(line);
						if (m.matches()) {
							matchcount++;
							accGoDb.put(m.group("accsynonym"), idGoDb.get(m.group("termid")));
						}
					}
					System.out.println(matchcount + " of " + linecount + " lines matched");
				}
				entry = tais.getNextTarEntry();
			}
		} finally {
			tais.close();
		}
		System.out.println("Finished reading and adding term synonyms");
	}
	
	// Check consistency of root term annotation count and information content
	private static void checkConsistency(){
		if (accGoDb.get(bpRootAcc).getFrequency() != bpCount) {
			System.out.println("Consistency error in GO frequencies of the biological process ontology");
			System.out.println("Term frequency: " + accGoDb.get(bpRootAcc).getFrequency());
		}
		if (accGoDb.get(ccRootAcc).getFrequency() != ccCount) {
			System.out.println("Consistency error in GO frequencies of the cellular component ontology");
			System.out.println("Term frequency: " + accGoDb.get(ccRootAcc).getFrequency());
		}
		if (accGoDb.get(mfRootAcc).getFrequency() != mfCount) {
			System.out.println("Consistency error in GO frequencies of the molecular function ontology");
			System.out.println("Term frequency: " + accGoDb.get(mfRootAcc).getFrequency());
		}
		if (accGoDb.get(bpRootAcc).getInformationContent() != 0.0) {
			System.out.println("Root term of biological process ontology has a non zero information content: " + accGoDb.get(bpRootAcc).getInformationContent());
		}
		if (accGoDb.get(ccRootAcc).getInformationContent() != 0.0) {
			System.out.println("Root term of cellular component has a non zero information content: " + accGoDb.get(ccRootAcc).getInformationContent());
		}
		if (accGoDb.get(mfRootAcc).getInformationContent() != 0.0) {
			System.out.println("Root term of molecular function ontology has a non zero information content: " + accGoDb.get(mfRootAcc).getInformationContent());
		}
	}
}
