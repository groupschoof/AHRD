package ahrd.model;

import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Utils.getAHRDdir;

import java.io.*;
import java.net.URL;
import java.nio.channels.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import javax.annotation.Nonnull;

import org.apache.commons.compress.archivers.tar.*;
import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.formats.FunctionalSyntaxDocumentFormat;
import org.semanticweb.owlapi.model.IRI;
import org.semanticweb.owlapi.model.OWLAnnotation;
import org.semanticweb.owlapi.model.OWLAnnotationAssertionAxiom;
import org.semanticweb.owlapi.model.OWLAnnotationProperty;
import org.semanticweb.owlapi.model.OWLAnnotationValue;
import org.semanticweb.owlapi.model.OWLClass;
import org.semanticweb.owlapi.model.OWLDataFactory;
import org.semanticweb.owlapi.model.OWLDocumentFormat;
import org.semanticweb.owlapi.model.OWLEntity;
import org.semanticweb.owlapi.model.OWLLiteral;
import org.semanticweb.owlapi.model.OWLNamedIndividual;
import org.semanticweb.owlapi.model.OWLOntology;
import org.semanticweb.owlapi.model.OWLOntologyCreationException;
import org.semanticweb.owlapi.model.OWLOntologyManager;
import org.semanticweb.owlapi.model.OWLOntologyStorageException;
import org.semanticweb.owlapi.reasoner.ConsoleProgressMonitor;
import org.semanticweb.owlapi.reasoner.Node;
import org.semanticweb.owlapi.reasoner.NodeSet;
import org.semanticweb.owlapi.reasoner.OWLReasoner;
import org.semanticweb.owlapi.reasoner.OWLReasonerConfiguration;
import org.semanticweb.owlapi.reasoner.OWLReasonerFactory;
import org.semanticweb.owlapi.reasoner.SimpleConfiguration;
import org.semanticweb.owlapi.reasoner.structural.StructuralReasonerFactory;
import org.semanticweb.owlapi.search.EntitySearcher;
import org.semanticweb.owlapi.util.DefaultPrefixManager;

import java.lang.Math;


public class GOdatabase {

	private static final String reviewedUniProtURL = "ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz";
	private static final String reviewedUniProtFileName = "uniprot_sprot.dat.gz";
	private static final Pattern reviewedUniProtGoAnnotationRegex = Pattern.compile("^DR\\s{3}GO;\\s+(?<goTerm>GO:\\d{7}).*");
	private static final String geneOntologyMonthlyOwlReleaseURL = "http://current.geneontology.org/ontology/go-basic.owl";
	private static final String geneOntologyMonthlyOwlReleaseFileName = "go-basic.owl";
	private static final String geneOntologyOwlPrefix = "http://purl.obolibrary.org/obo/GO_";
	private static final Pattern geneOntologyMYSQLdumpTermTableRegex = Pattern.compile("^(?<id>\\d+)\\t(?<name>[^\\t]+)\\t(?<termtype>[^\\t]+)\\t(?<acc>[^\\t]+)\\t(?<isobsolete>\\d)\\t(?<isroot>\\d)\\t(?<isrelation>\\d)$");
	private static final Pattern geneOntologyMYSQLdumpGraphPathRegex = Pattern.compile("^(?<id>\\d+)\\t(?<term1id>\\d+)\\t(?<term2id>\\d+)\\t(?<relationshiptypeid>1|25)\\t(?<distance>\\d+)\\t(?<relationdistance>\\d+)$");
	private static final Pattern geneOntologyMYSQLdumpTermSynonymRegex = Pattern.compile("^(?<termid>\\d+)\\t(?<termsynonym>[^\\t]+)\\t(?<accsynonym>GO:\\d{7})\\t(?<synonymtypeid>[^\\t]+)\\t(?<synonymcategoryid>[^\\t]+)$");
	private static final String bpRootAcc = "GO:0008150";
	private static final String ccRootAcc = "GO:0005575";
	private static final String mfRootAcc = "GO:0003674";
	private static final String serializedAccGoDBFileName = "accGoDb.ser";
	private Map<String, GOterm> goDb = new HashMap<String, GOterm>();
	
	public GOdatabase() throws Exception {
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
		
		// New way to build the database using OWLAPI (https://github.com/owlcs/owlapi)
		Map<String, GOterm> goDbBuildFromOWL = buildGoDbFromOWL(pathToGoDatabase);
	}
	
	private HashMap<String, GOterm> buildGoDbFromOWL(String pathToGoDatabase) throws Exception {
		System.out.println("Building GO database using OWLAPI (https://github.com/owlcs/owlapi):");
		HashMap<String, GOterm> accGoDb = new HashMap<String, GOterm>();

		String reviewedUniProtFilePath =  pathToGoDatabase + reviewedUniProtFileName;		
		// Download SwissProt if not already on drive
		if (!new File(reviewedUniProtFilePath).exists()) {
			System.out.println("Downloading reviewed Uniprot (aprox. 550MB) from:\n"+ reviewedUniProtURL); 
			download(reviewedUniProtURL, reviewedUniProtFilePath);
		}

		String geneOntologyMonthlyOwlReleasePath = pathToGoDatabase + geneOntologyMonthlyOwlReleaseFileName;
		// Download gene ontology OWL if not alrady on drive
		if (!new File(geneOntologyMonthlyOwlReleasePath).exists()) {
			System.out.println("Downloading GO database (aprox. 150MB) from: " + geneOntologyMonthlyOwlReleaseURL);
			download(geneOntologyMonthlyOwlReleaseURL, geneOntologyMonthlyOwlReleasePath);
		}
				
		// Get hold of a manager to work with
		OWLOntologyManager manager = OWLManager.createOWLOntologyManager();
		// Load the gene ontology from file
		OWLOntology ontology = manager.loadOntologyFromOntologyDocument(new File(geneOntologyMonthlyOwlReleasePath));
		// Print the ontology ID
		System.out.println("Loaded ontology: " + ontology.getOntologyID());
		// Print the format of the ontology
		System.out.println("Ontology format:" + manager.getOntologyFormat(ontology).getKey());
		// An OWLReasoner provides the basic query functionality needed, for example the ability
        // to obtain the subclasses of a class etc. To do this reasoner factory is used.
		OWLReasonerFactory reasonerFactory = new StructuralReasonerFactory();
		OWLReasoner reasoner = reasonerFactory.createReasoner(ontology);
		// Ask the reasoner to do all the necessary work now
        reasoner.precomputeInferences();
        // Determine if the ontology is actually consistent.
        System.out.println("Consistent: " + reasoner.isConsistent());
		//OWLDataFactory dataFactory = manager.getOWLDataFactory();

		// Iterate over all classes and store terms in map
		for (OWLClass cls : ontology.getClassesInSignature()) {
			String id = "";
			String label = "";
			String nameSpace = "";
	        for(OWLAnnotation axiom : EntitySearcher.getAnnotations(cls, ontology)) {
	    		if(axiom.getProperty().getIRI().equals(IRI.create("http://www.geneontology.org/formats/oboInOwl#id"))) {
	    			if(axiom.getValue() instanceof OWLLiteral) {
	    	            OWLLiteral val = (OWLLiteral) axiom.getValue();
	    	            id = val.getLiteral();
	    	        }
	    		}
	    		if(axiom.getProperty().getIRI().equals(IRI.create("http://www.w3.org/2000/01/rdf-schema#label"))) {
	    			if(axiom.getValue() instanceof OWLLiteral) {
	    	            OWLLiteral val = (OWLLiteral) axiom.getValue();
	    	            label = val.getLiteral();
	    	        }
	    		}
	    		if(axiom.getProperty().getIRI().equals(IRI.create("http://www.geneontology.org/formats/oboInOwl#hasOBONamespace"))) {
	    			if(axiom.getValue() instanceof OWLLiteral) {
	    	            OWLLiteral val = (OWLLiteral) axiom.getValue();
	    	            nameSpace = val.getLiteral();
	    	        }
	    		}
	    	}
	        if(!id.equals("")) {
	        	accGoDb.put(id, new GOterm(id, label, nameSpace));
	        }
        }
		System.out.println("Size of accGoDb based on OWL: " + accGoDb.size());
		// Iterate over all classes again and store ancestry
		Integer ancestryCounter = 0;
		for (OWLClass cls : ontology.getClassesInSignature()) {
			String childId = "";
			for(OWLAnnotation axiom : EntitySearcher.getAnnotations(cls, ontology)) {
	    		if(axiom.getProperty().getIRI().equals(IRI.create("http://www.geneontology.org/formats/oboInOwl#id"))) {
	    			if(axiom.getValue() instanceof OWLLiteral) {
	    	            OWLLiteral val = (OWLLiteral) axiom.getValue();
	    	            childId = val.getLiteral();
	    	        }
	    		}
			}
			if (!childId.equals("")) {
				Set<GOterm> ancestry = new HashSet<GOterm>();
		    	Set<OWLClass> superClses = reasoner.getSuperClasses(cls, false).getFlattened();
		    	for (OWLClass superCls : superClses) {
		    		for(OWLAnnotation axiom : EntitySearcher.getAnnotations(superCls, ontology)) {
			    		if(axiom.getProperty().getIRI().equals(IRI.create("http://www.geneontology.org/formats/oboInOwl#id"))) {
			    			if(axiom.getValue() instanceof OWLLiteral) {
			    	            OWLLiteral val = (OWLLiteral) axiom.getValue();
			    	            ancestry.add(accGoDb.get(val.getLiteral()));
			    	        }
			    		}
		    		}
		    	}
		    	ancestryCounter += ancestry.size();
		    	accGoDb.get(childId).setAncestry(ancestry);
			}
		}
		
		System.out.println("Number of ancestry terms: " + ancestryCounter);
		return accGoDb;
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

	private HashMap<String, GOterm> buildGoDbFromFile(String pathToGoDatabase) throws Exception {
		System.out.println("Building GO database:");
		HashMap<String, GOterm> accGoDb = new HashMap<String, GOterm>();

		String reviewedUniProtFilePath =  pathToGoDatabase + reviewedUniProtFileName;		
		// Download SwissProt if not already on drive
		if (!new File(reviewedUniProtFilePath).exists()) {
			System.out.println("Downloading reviewed Uniprot (aprox. 550MB) from:\n"+ reviewedUniProtURL); 
			download(reviewedUniProtURL, reviewedUniProtFilePath);
		}

		String geneOntologyMonthlyOwlReleasePath = pathToGoDatabase + geneOntologyMonthlyOwlReleaseFileName;
		// Download gene ontology OWL if not alrady on drive
		if (!new File(geneOntologyMonthlyOwlReleasePath).exists()) {
			System.out.println("Downloading GO database (aprox. 150MB) from:\n" + geneOntologyMonthlyOwlReleaseURL);
			download(geneOntologyMonthlyOwlReleaseURL, geneOntologyMonthlyOwlReleasePath);
		}		
		
		Map<Integer, GOterm> idGoDb = new HashMap<Integer, GOterm>();
		
		//////////////////////////////////////////////////////////////////////////////////////
		TarArchiveInputStream tais = null;
		// Read term table and fill ID based GO term database
		try {
			tais = new TarArchiveInputStream(new GZIPInputStream(new FileInputStream(geneOntologyMonthlyOwlReleasePath)));
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
			tais = new TarArchiveInputStream(new GZIPInputStream(new FileInputStream(geneOntologyMonthlyOwlReleasePath)));
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
			tais = new TarArchiveInputStream(new GZIPInputStream(new FileInputStream(geneOntologyMonthlyOwlReleasePath)));
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
									parent.setFrequency(parent.getFrequency() + 1);
									}
								break;
							case "cellular_component":
								ccCount++;
								for(GOterm parent : ancestry) {
									parent.setFrequency(parent.getFrequency() + 1);
								}
								break;
							case "molecular_function":
								mfCount++;
								for(GOterm parent : ancestry) {
									parent.setFrequency(parent.getFrequency() + 1);
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
	
		// Calculate the information content from the term probabilities
		mapIter = idGoDb.entrySet().iterator();
		term = null;
		int metacount = 0;
		while (mapIter.hasNext()) {
			Map.Entry<Integer, GOterm> entry = mapIter.next();
			term = entry.getValue();
			switch(term.getOntology()) {
			case "biological_process":
				term.setProbability((double)term.getFrequency()/bpCount);
				term.setInformationContent(-1*Math.log(term.getProbability()));
				break;
			case "cellular_component":
				term.setProbability((double)term.getFrequency()/ccCount);
				term.setInformationContent(-1*Math.log(term.getProbability()));
				break;
			case "molecular_function":
				term.setProbability((double)term.getFrequency()/mfCount);
				term.setInformationContent(-1*Math.log(term.getProbability()));
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
