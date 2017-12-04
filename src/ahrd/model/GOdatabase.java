package ahrd.model;

import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Utils.getAHRDdir;

import java.io.*;
import java.net.URL;
import java.nio.channels.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import org.apache.commons.compress.archivers.tar.*;
import java.lang.Math;

public class GOdatabase {

	private static final String reviewedUniProtURL = "ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz";
	private static final String reviewedUniProtFileName = "uniprot_sprot.dat.gz";
	private static final Pattern reviewedUniProtGoAnnotationRegex = Pattern.compile("^DR\\s{3}GO;\\s+(?<goTerm>GO:\\d{7}).*");
	private static final String geneOntologyMYSQLdumpURL = "http://archive.geneontology.org/latest-termdb/go_daily-termdb-tables.tar.gz";
	private static final String geneOntologyMYSQLdumpFileName = "go_daily-termdb-tables.tar.gz";
	private static final Pattern geneOntologyMYSQLdumpTermTableRegex = Pattern.compile("^(?<id>\\d+)\\t(?<name>[^\\t]+)\\t(?<termtype>[^\\t]+)\\t(?<acc>[^\\t]+)\\t(?<isobsolete>\\d)\\t(?<isroot>\\d)\\t(?<isrelation>\\d)$");
	private static final Pattern geneOntologyMYSQLdumpGraphPathRegex = Pattern.compile("^(?<id>\\d+)\\t(?<term1id>\\d+)\\t(?<term2id>\\d+)\\t(?<relationshiptypeid>1|25)\\t(?<distance>\\d+)\\t(?<relationdistance>\\d+)$");
	private static final Pattern geneOntologyMYSQLdumpTermSynonymRegex = Pattern.compile("^(?<termid>\\d+)\\t(?<termsynonym>[^\\t]+)\\t(?<accsynonym>GO:\\d{7})\\t(?<synonymtypeid>[^\\t]+)\\t(?<synonymcategoryid>[^\\t]+)$");
	private static final String bpRootAcc = "GO:0008150";
	private static final String ccRootAcc = "GO:0005575";
	private static final String mfRootAcc = "GO:0003674";
	private static final String serializedAccGoDBFileName = "accGoDb.ser";
	private Map<String, GOterm> goDb = new HashMap<String, GOterm>();
	
	public GOdatabase() throws FileNotFoundException, IOException {
		String pathToGoDatabase = new String();
		if (getSettings().getPathToGoDatabase() != null) {
			pathToGoDatabase = getSettings().getPathToGoDatabase();
		} else {
			pathToGoDatabase = getAHRDdir(this.getClass()) + "data/";
		}
		if (!pathToGoDatabase.endsWith("/")) {
			pathToGoDatabase = pathToGoDatabase + "/";
		}
		String serializedAccGoDBFilePath = pathToGoDatabase + serializedAccGoDBFileName;
		if (new File(serializedAccGoDBFilePath).exists()) {
			goDb = deserializeAccGoDb(serializedAccGoDBFilePath);
		} else {
			goDb = buildGoDbFromFile(pathToGoDatabase);
		}
	}
	
	// Deserialize previously generated GO database
	@SuppressWarnings("unchecked")
	private HashMap<String, GOterm> deserializeAccGoDb(String serializedAccGoDBFilePath) {
		if (new File(serializedAccGoDBFilePath).exists()) {
			try {
				FileInputStream fileIn = new FileInputStream(serializedAccGoDBFilePath);
				ObjectInputStream in = new ObjectInputStream(fileIn);
				HashMap<String, GOterm> goDbFromFile = (HashMap<String, GOterm>) in.readObject();
				in.close();
				fileIn.close();
				System.out.println("GO database sucessfully deserialized. Size: " + goDbFromFile.size());
				return goDbFromFile;
			} catch (java.io.InvalidClassException exp) {
				System.out.println("GOterm class incompatibility discovered while deserializing GO database object."
						+ "\nThe accGoDb.ser file in the ./AHRD/data folder seems to have been created with a differnt build of AHRD."
						+ "\nPlease trigger the automatic rebuild process for the GO database by deleting the accGoDb.ser file before you run AHRD again.");
				exp.printStackTrace();
				return null;
			} catch (IOException i) {
				i.printStackTrace();
				return null;
			} catch (ClassNotFoundException c) {
				System.out.println("GOterm class not found for deserialization of GO database object");
				c.printStackTrace();
				return null;
			}
		}
		return null;
	}

	private HashMap<String, GOterm> buildGoDbFromFile(String pathToGoDatabase) throws FileNotFoundException, IOException {
		System.out.println("Building GO database:");
		HashMap<String, GOterm> accGoDb = new HashMap<String, GOterm>();

		String reviewedUniProtFilePath =  pathToGoDatabase + reviewedUniProtFileName;		
		// Download SwissProt if not already on drive
		if (!new File(reviewedUniProtFilePath).exists()) {
			System.out.println("Downloading reviewed Uniprot (aprox. 550MB) from:\n"+ reviewedUniProtURL); 
			download(reviewedUniProtURL, reviewedUniProtFilePath);
		}

		String geneOntologyMYSQLdumpFilePath = pathToGoDatabase + geneOntologyMYSQLdumpFileName;
		// Download gene ontology mysql data base dump if not alrady on drive
		if (!new File(geneOntologyMYSQLdumpFilePath).exists()) {
			System.out.println("Downloading GO database (aprox. 12MB) from:\n" + geneOntologyMYSQLdumpURL);
			download(geneOntologyMYSQLdumpURL, geneOntologyMYSQLdumpFilePath);
		}
		
		Map<Integer, GOterm> idGoDb = new HashMap<Integer, GOterm>();
		
		TarArchiveInputStream tais = null;
		// Read term table and fill ID based GO term database
		try {
			tais = new TarArchiveInputStream(new GZIPInputStream(new FileInputStream(geneOntologyMYSQLdumpFilePath)));
			TarArchiveEntry entry = tais.getNextTarEntry();
			while (entry != null) {
				if (entry.getName().equals("go_daily-termdb-tables/term.txt")) {
					BufferedReader termtablebr = null;
					termtablebr = new BufferedReader(new InputStreamReader(tais));
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
				entry = tais.getNextTarEntry();
			}
		} finally {
			tais.close();
		}

		// Read graph path table and fill the terms ancestry in goDB
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
							GOterm term = idGoDb.get(Integer.parseInt(m.group("term2id")));
							GOterm parent = idGoDb.get(Integer.parseInt(m.group("term1id")));
							if (term != null && parent != null) {
								term.addTermToAncestry(parent);
							}
						}
					}
					System.out.println(matchcount + " of " + linecount + " graph path lines matched.");
				}
				entry = tais.getNextTarEntry();
			}
		} finally {
			tais.close();
		}

		// Build an accession based GO database
		Iterator<Map.Entry<Integer, GOterm>> mapIter = idGoDb.entrySet().iterator();
		while (mapIter.hasNext()) {
			Map.Entry<Integer, GOterm> entry = mapIter.next();
			accGoDb.put(entry.getValue().getAccession(), entry.getValue());
		}
	
		// Read in term synonyms and add them to accession based GO database
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
							accGoDb.put(m.group("accsynonym"), idGoDb.get(Integer.parseInt(m.group("termid"))));
						}
					}
					System.out.println(matchcount + " of " + linecount + " synonym table lines matched.");
				}
				entry = tais.getNextTarEntry();
			}
		} finally {
			tais.close();
		}

		// Read in SwissProt annotations
		// Counts each annotation towards a GO term and its complete ancestry  
		System.out.println("Reading and counting aprox 2.7 million SwissProt annotations ...");
		int bpCount = 0;
		int ccCount = 0;
		int mfCount = 0;
		String line, termAcc= new String();
		Integer linecount = 0;
		Integer matchcount = 0;
		BufferedReader br = null;
		GOterm term = null;
		Set<GOterm> ancestry = null;
		try {
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(reviewedUniProtFilePath))));
			while ((line = br.readLine()) != null) {
				linecount++;
				Matcher m = reviewedUniProtGoAnnotationRegex.matcher(line);
				if (m.matches()) {
					matchcount++;
					termAcc = m.group("goTerm");
					term = accGoDb.get(termAcc);
					if (term == null) {
						System.out.println("Error: GO-accession of SwissProt annotation not found in Gene Ontology.");
					} else {
						ancestry = term.getAncestry();
						if (!term.getObsolete()) {
							switch(term.getOntology()) {
							case "biological_process":
								bpCount++;
								for(GOterm parent : ancestry) {
									parent.setAnnotationCount(parent.getAnnotationCount() + 1);
									}
								break;
							case "cellular_component":
								ccCount++;
								for(GOterm parent : ancestry) {
									parent.setAnnotationCount(parent.getAnnotationCount() + 1);
								}
								break;
							case "molecular_function":
								mfCount++;
								for(GOterm parent : ancestry) {
									parent.setAnnotationCount(parent.getAnnotationCount() + 1);
								}
								break;
							default:
								break; //meta term: Not counted towards any of the three ontologies
							}

						}
					}
				}
			}
		} finally {
			br.close();
		}
		System.out.println(matchcount + " of " + linecount + " SwissProt lines matched.");
		System.out.println("bpCount: " + bpCount);
		System.out.println("ccCount: " + ccCount);
		System.out.println("mfCount: " + mfCount);
	
		// Calculate the information content from the term annotation frequencies
		mapIter = idGoDb.entrySet().iterator();
		term = null;
		int metacount = 0;
		while (mapIter.hasNext()) {
			Map.Entry<Integer, GOterm> entry = mapIter.next();
			term = entry.getValue();
			switch(term.getOntology()) {
			case "biological_process":
				term.setAnnotationFrequency((double)term.getAnnotationCount()/bpCount);
				term.setInformationContent(-1*Math.log(term.getAnnotationFrequency()));
				break;
			case "cellular_component":
				term.setAnnotationFrequency((double)term.getAnnotationCount()/ccCount);
				term.setInformationContent(-1*Math.log(term.getAnnotationFrequency()));
				break;
			case "molecular_function":
				term.setAnnotationFrequency((double)term.getAnnotationCount()/mfCount);
				term.setInformationContent(-1*Math.log(term.getAnnotationFrequency()));
				break;
			default:
				metacount++;//meta term: Not counted towards any of the three ontologies; omit storing in accession based term database  
				break;
			}
			
		}
		System.out.println("metacount: " + metacount);
		System.out.println("Finished calculating information content.");
		System.out.println("Size of accGoDb: " + accGoDb.size());
		
		// Check consistency of root term annotation count and information content
		if (accGoDb.get(bpRootAcc).getAnnotationCount() != bpCount) {
			System.out.println("Consistency error in GO frequencies of the biological process ontology");
			System.out.println("Term frequency: " + accGoDb.get(bpRootAcc).getAnnotationCount());
		}
		if (accGoDb.get(ccRootAcc).getAnnotationCount() != ccCount) {
			System.out.println("Consistency error in GO frequencies of the cellular component ontology");
			System.out.println("Term frequency: " + accGoDb.get(ccRootAcc).getAnnotationCount());
		}
		if (accGoDb.get(mfRootAcc).getAnnotationCount() != mfCount) {
			System.out.println("Consistency error in GO frequencies of the molecular function ontology");
			System.out.println("Term frequency: " + accGoDb.get(mfRootAcc).getAnnotationCount());
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

		// Serialize accession based GO database and write it to file	
		try {
			String serializedAccGoDBFilePath = pathToGoDatabase + serializedAccGoDBFileName;
			FileOutputStream fileOut = new FileOutputStream(serializedAccGoDBFilePath);
			ObjectOutputStream out = new ObjectOutputStream(fileOut);
			out.writeObject(accGoDb);
			out.close();
			fileOut.close();
			System.out.println("Serialized data (aprox. 11MB) saved in: " + serializedAccGoDBFilePath);
		} catch (IOException i) {
			i.printStackTrace();
		}
		System.out.println("Done building new GO database!");
		return accGoDb;
	}
	
	private static void download(String url, String file) throws IOException {
		URL website = new URL(url);
		ReadableByteChannel rbc = Channels.newChannel(website.openStream());
		FileOutputStream fos = new FileOutputStream(file);
		fos.getChannel().transferFrom(rbc, 0, Long.MAX_VALUE);
		fos.close();
	}

	public static String getBprootacc() {
		return bpRootAcc;
	}

	public static String getCcrootacc() {
		return ccRootAcc;
	}

	public static String getMfrootacc() {
		return mfRootAcc;
	}

	public Map<String, GOterm> getMap() {
		return goDb;
	}

}
