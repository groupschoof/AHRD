package ahrd.model;

import static ahrd.controller.Settings.getSettings;
import static ahrd.controller.Utils.getAHRDdir;
import static ahrd.controller.Utils.roundToNDecimalPlaces;

import java.io.*;
import java.net.URL;
import java.nio.channels.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.semanticweb.owlapi.apibinding.OWLManager;
import org.semanticweb.owlapi.model.IRI;
import org.semanticweb.owlapi.model.OWLAnnotation;
import org.semanticweb.owlapi.model.OWLClass;
import org.semanticweb.owlapi.model.OWLDataFactory;
import org.semanticweb.owlapi.model.OWLLiteral;
import org.semanticweb.owlapi.model.OWLOntology;
import org.semanticweb.owlapi.model.OWLOntologyCreationException;
import org.semanticweb.owlapi.model.OWLOntologyManager;
import org.semanticweb.owlapi.reasoner.OWLReasoner;
import org.semanticweb.owlapi.reasoner.OWLReasonerFactory;
import org.semanticweb.owlapi.reasoner.structural.StructuralReasonerFactory;
import org.semanticweb.owlapi.search.EntitySearcher;

import java.lang.Math;


public class GOdatabase {

	private static final String reviewedUniProtURL = "ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz";
	private static final String reviewedUniProtFileName = "uniprot_sprot.dat.gz";
	private static final Pattern reviewedUniProtGoAnnotationRegex = Pattern.compile("^DR\\s{3}GO;\\s+(?<goTerm>GO:\\d{7}).*");
	private static final String geneOntologyMonthlyOwlReleaseURL = "http://current.geneontology.org/ontology/go-basic.owl";
	private static final String geneOntologyMonthlyOwlReleaseFileName = "go-basic.owl";
	private static final String bpRootAcc = "GO:0008150";
	private static final String ccRootAcc = "GO:0005575";
	private static final String mfRootAcc = "GO:0003674";

	private static final String serializedAccGoDbOwlFileName = "accGoDbOwl.ser";
	private Map<String, GOterm> goDb = new HashMap<String, GOterm>();
	
	public GOdatabase() throws OWLOntologyCreationException, IOException {
		String pathToGoDatabase = new String();
		if (getSettings().getPathToGoDatabase() != null) {
			pathToGoDatabase = getSettings().getPathToGoDatabase();
		} else {
			pathToGoDatabase = getAHRDdir(this.getClass()) + "data/";
		}
		if (!pathToGoDatabase.endsWith("/")) {
			pathToGoDatabase = pathToGoDatabase + "/";
		}
		// New way to build the database using OWLAPI (https://github.com/owlcs/owlapi)
		String serializedAccGoDbOwlFilePath = pathToGoDatabase + serializedAccGoDbOwlFileName;
		if (new File(serializedAccGoDbOwlFilePath).exists()) {
			goDb = deserializeAccGoDb(serializedAccGoDbOwlFilePath);
		} else {
			goDb = buildGoDbFromOWL(pathToGoDatabase);
		}
	}
	
	private HashMap<String, GOterm> buildGoDbFromOWL(String pathToGoDatabase) throws IOException, OWLOntologyCreationException  {
		System.out.println("Building GO database using OWLAPI (https://github.com/owlcs/owlapi)");
		HashMap<String, GOterm> accGoDb = new HashMap<String, GOterm>();

		String reviewedUniProtFilePath =  pathToGoDatabase + reviewedUniProtFileName;		
		// Download SwissProt if not already on drive
		if (!new File(reviewedUniProtFilePath).exists()) {
			System.out.println("Downloading reviewed Uniprot (aprox. 570MB) from:\n"+ reviewedUniProtURL); 
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
		System.out.println("Loaded ontology version: " + ontology.getOntologyID().getVersionIRI().get());
		// Print the format of the ontology
		System.out.println("Ontology format: " + manager.getOntologyFormat(ontology).getKey());
		// An OWLReasoner provides the basic query functionality needed, for example the ability
        // to obtain the subclasses of a class etc. To do this reasoner factory is used.
		OWLReasonerFactory reasonerFactory = new StructuralReasonerFactory();
		OWLReasoner reasoner = reasonerFactory.createReasoner(ontology);
		// Ask the reasoner to do all the necessary work now
        reasoner.precomputeInferences();
        // Determine if the ontology is actually consistent.
        System.out.println("Consistent: " + reasoner.isConsistent());
		OWLDataFactory dataFactory = manager.getOWLDataFactory();

		// Iterate over all classes and store terms in map
		long processGoStartTime = System.currentTimeMillis();
		for (OWLClass cls : ontology.getClassesInSignature()) {
			String id = "";
			Set<String> altIds = new HashSet<String>();
			String label = "";
			String nameSpace = "";
			Boolean deprecated = false;
	        for(OWLAnnotation annotation : EntitySearcher.getAnnotations(cls, ontology)) {
	        	String annotationIri = annotation.getProperty().getIRI().toString();
	        	if(annotationIri.equals("http://www.geneontology.org/formats/oboInOwl#id")) {
	    			if(annotation.getValue() instanceof OWLLiteral) {
	    	            OWLLiteral val = (OWLLiteral) annotation.getValue();
	    	            id = val.getLiteral();
	    	        }
	    		}
	    		if(annotationIri.equals("http://www.geneontology.org/formats/oboInOwl#hasAlternativeId")) {
	    			if(annotation.getValue() instanceof OWLLiteral) {
	    	            OWLLiteral val = (OWLLiteral) annotation.getValue();
	    	            altIds.add(val.getLiteral());
	    	        }
	    		}
	    		if(annotationIri.equals("http://www.w3.org/2000/01/rdf-schema#label")) {
	    			if(annotation.getValue() instanceof OWLLiteral) {
	    	            OWLLiteral val = (OWLLiteral) annotation.getValue();
	    	            label = val.getLiteral();
	    	        }
	    		}
	    		if(annotationIri.equals("http://www.geneontology.org/formats/oboInOwl#hasOBONamespace")) {
	    			if(annotation.getValue() instanceof OWLLiteral) {
	    	            OWLLiteral val = (OWLLiteral) annotation.getValue();
	    	            nameSpace = val.getLiteral();
	    	        }
	    		}
	        	if(annotationIri.equals("http://www.w3.org/2002/07/owl#deprecated")) {
	        		if(annotation.getValue() instanceof OWLLiteral) {
	        			OWLLiteral val = (OWLLiteral) annotation.getValue();
	        			deprecated = val.parseBoolean();
	        		}
	        	}
	    	}
	        if (!id.equals("")) {
	        	GOterm term = new GOterm(id, label, nameSpace, deprecated);
	        	accGoDb.put(id, term);
	        	for (String altId : altIds) {
	        		accGoDb.put(altId, term);
	        	}
	        }
        }
		// Iterate over all classes again and store ancestry
		Integer ancestryCounter = 0;
		for (OWLClass cls : ontology.getClassesInSignature()) {
			String id = "";
			for(OWLAnnotation annotation : EntitySearcher.getAnnotations(cls, ontology, dataFactory.getOWLAnnotationProperty(IRI.create("http://www.geneontology.org/formats/oboInOwl#id")))) {
	    			if(annotation.getValue() instanceof OWLLiteral) {
	    	            OWLLiteral val = (OWLLiteral) annotation.getValue();
	    	            id = val.getLiteral();
	    	        }
			}
			if (!id.equals("")){
				Set<GOterm> ancestry = new HashSet<GOterm>();
		    	Set<OWLClass> superClses = reasoner.getSuperClasses(cls, false).getFlattened();
		    	for (OWLClass superCls : superClses) {
		    		for(OWLAnnotation annotation : EntitySearcher.getAnnotations(superCls, ontology, dataFactory.getOWLAnnotationProperty(IRI.create("http://www.geneontology.org/formats/oboInOwl#id")))) {
		    			if(annotation.getValue() instanceof OWLLiteral) {
		    	            OWLLiteral val = (OWLLiteral) annotation.getValue();
		    	            ancestry.add(accGoDb.get(val.getLiteral()));
		    	        }
		    		}
		    	}
	    	   	if (accGoDb.get(id) == null) {
		    		System.err.println("Id " + id + " not found in accGoDb!");
		    	} else {
		    		/* Add the go term itself to its own ancestry.
		    		 * It was this way in the deprecated mysql dumps and is kept up this way here, to not break other methods using the GOdatabase.
		    		 * (a particular example: ahrd.model.EvaluationScoreCalculator.calcSemSimGoAnnotationScore())  
		    		 */
		    		ancestry.add(accGoDb.get(id));
		    		ancestryCounter += ancestry.size();
		    		accGoDb.get(id).setAncestry(ancestry);
		    	}
	    	}
		}
		long processGoFinishTime = System.currentTimeMillis();
		long processGoTimeElapsedMilli = processGoFinishTime - processGoStartTime;
		double processGoTimeElapsedSec = roundToNDecimalPlaces((double) processGoTimeElapsedMilli / (double) 1000, 1);

		System.out.println("Number of ancestry relations retrieved from OWL: " + ancestryCounter);
		System.out.println("Time to process the ancestry: " + processGoTimeElapsedSec + " s");
		
		// Read in SwissProt annotations
		// Counts each annotation towards a GO term and its complete ancestry  
		System.out.println("Reading and counting aprox. 3 million SwissProt annotations ...");
		long processSwissProtStartTime = System.currentTimeMillis();
		int bpCount = 0;
		int ccCount = 0;
		int mfCount = 0;
		String line = new String();
		Integer linecount = 0;
		Integer matchcount = 0;
		BufferedReader br = null;
		try {
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(reviewedUniProtFilePath))));
			while ((line = br.readLine()) != null) {
				linecount++;
				Matcher m = reviewedUniProtGoAnnotationRegex.matcher(line);
				if (m.matches()) {
					matchcount++;
					String termAcc = m.group("goTerm");
					GOterm term = accGoDb.get(termAcc);
					if (term == null) {
						System.err.println("Error: GO-accession (" + termAcc + ") of SwissProt annotation not found in Gene Ontology.");
					} else {
						Set<GOterm> ancestry = term.getAncestry();
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
		System.out.println("Found " + matchcount + " go annotation lines in " + linecount + " lines of SwissProt dat file.");
		System.out.println("Biological process terms: " + bpCount);
		System.out.println("Cellular component terms: " + ccCount);
		System.out.println("Molecular function terms: " + mfCount);
		long processSwissProtFinishTime = System.currentTimeMillis();
		long processSwissProtElapsedTimeMilli = processSwissProtFinishTime - processSwissProtStartTime;
		double proccesSwissProtElapsedTimeSec = roundToNDecimalPlaces((double) processSwissProtElapsedTimeMilli /(double) 1000, 1);
		System.out.println("Time to process SwissProt annotations: " + proccesSwissProtElapsedTimeSec + " s");
		// Calculate the information content from the term annotation frequencies
		for (GOterm term : accGoDb.values() ){
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
				break; //meta term: Not counted towards any of the three ontologies
			}
		}
		System.out.println("Finished calculating information content.");
		System.out.println("GO terms in database: " + accGoDb.size());
		
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
			System.err.println("Root term of biological process ontology has a non zero information content: " + accGoDb.get(bpRootAcc).getInformationContent());
		}
		if (accGoDb.get(ccRootAcc).getInformationContent() != 0.0) {
			System.err.println("Root term of cellular component has a non zero information content: " + accGoDb.get(ccRootAcc).getInformationContent());
		}
		if (accGoDb.get(mfRootAcc).getInformationContent() != 0.0) {
			System.err.println("Root term of molecular function ontology has a non zero information content: " + accGoDb.get(mfRootAcc).getInformationContent());
		}

		// Serialize accession based GO database and write it to file	
		try {
			String serializedAccGoDBFilePath = pathToGoDatabase + serializedAccGoDbOwlFileName;
			FileOutputStream fileOut = new FileOutputStream(serializedAccGoDBFilePath);
			ObjectOutputStream out = new ObjectOutputStream(fileOut);
			out.writeObject(accGoDb);
			out.close();
			fileOut.close();
			System.out.println("Serialized database (aprox. 10MB) saved in: " + serializedAccGoDBFilePath);
		} catch (IOException i) {
			i.printStackTrace();
		}
		System.out.println("Done building new GO database based on OWL file!");
		Double maxFiniteInfoContent = 0.0;
		GOterm maxFiniteInfoContentTerm = null;
		Double minInfoContent = Double.MAX_VALUE;
		GOterm minInfoContentTerm = null;
		for (GOterm term : accGoDb.values()) {
			if (Double.isFinite(term.getInformationContent()) && term.getInformationContent() > maxFiniteInfoContent) {
				maxFiniteInfoContent = term.getInformationContent();
				maxFiniteInfoContentTerm = term;
			}
			if (term.getInformationContent() > 0 && term.getInformationContent() < minInfoContent) {
				minInfoContent = term.getInformationContent();
				minInfoContentTerm = term;
			}
		}
		System.out.println("MinInfoContent:" + minInfoContentTerm.getAccession() + " " + minInfoContentTerm.getOntology() + " " + minInfoContentTerm.getAnnotationCount() + " " + minInfoContent);
		System.out.println("MaxInfoContent:" + maxFiniteInfoContentTerm.getAccession() + " " + maxFiniteInfoContentTerm.getOntology() + " " + maxFiniteInfoContentTerm.getAnnotationCount() + " " + maxFiniteInfoContent);
		int singletsCCOcount = 0;
		int singletsMFOcount = 0;
		int singletsBPOcount = 0;
		for (GOterm term : accGoDb.values()) {
			if (term.getAnnotationCount() == 1) {
				switch(term.getOntology()) {
				case "biological_process":
					singletsBPOcount++;
					break;
				case "cellular_component":
					singletsCCOcount++;
					break;
				case "molecular_function":
					singletsMFOcount++;
					break;
					}
			}
		}
		System.out.println("Single annotation terms counts:");
		System.out.println("BPO: " + singletsBPOcount);
		System.out.println("CCO: " + singletsCCOcount);
		System.out.println("MFO: " + singletsMFOcount);
		return accGoDb;
	}

	// Deserialize previously generated GO database
	@SuppressWarnings("unchecked")
	private HashMap<String, GOterm> deserializeAccGoDb(String serializedAccGoDBFilePath) {
		if (new File(serializedAccGoDBFilePath).exists()) {
			System.out.println("Proviously build GO database found: " + serializedAccGoDBFilePath);
			try {
				FileInputStream fileIn = new FileInputStream(serializedAccGoDBFilePath);
				ObjectInputStream in = new ObjectInputStream(fileIn);
				HashMap<String, GOterm> goDbFromFile = (HashMap<String, GOterm>) in.readObject();
				in.close();
				fileIn.close();
				System.out.println("Sucessfully deserialized. Size: " + goDbFromFile.size() + " terms");
				return goDbFromFile;
			} catch (java.io.InvalidClassException exp) {
				System.out.println("GOterm class incompatibility discovered while deserializing GO database object."
						+ "\nThe accGoDb.ser file in the ./AHRD/data folder seems to have been created with a differnt build of AHRD."
						+ "\nPlease trigger the automatic rebuild process for the GO database by deleting the accGoDbOwl.ser file before you run AHRD again.");
				exp.printStackTrace();
				return null;
			} catch (IOException i) {
				i.printStackTrace();
				return null;
			} catch (ClassNotFoundException c) {
				System.err.println("GOterm class not found for deserialization of GO database object");
				c.printStackTrace();
				return null;
			}
		}
		return null;
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
