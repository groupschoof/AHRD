package ahrd.model;

import static ahrd.controller.Settings.SHORT_ACCESSION_GROUP_NAME;
import static ahrd.controller.Settings.GO_TERM_GROUP_NAME;
import static ahrd.controller.Settings.getSettings;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GoAnnotationReference {
	/**
	 * Parses the tabular reference Gene Ontology term annotations (GOA) for
	 * proteins in the searched Blast databases. These GOAs will then be used to
	 * annotate the query proteins with GO terms. The important restriction is,
	 * that only those reference GO annotations will be extracted that match one
	 * of the BlastResults found in the respective Blast searches
	 * <code>uniqueShortAccessions</code>.
	 * 
	 * @param Set
	 *            <String> uniqueShortAccessions - The unique short accessions
	 *            of BlastResults found in the respective Blast searches
	 * @return Map<String, Set<String>> - BlastResult short-accessions mapped to
	 *         Sets of GO terms
	 * @throws IOException
	 */
	public static Map<String, Set<String>> parseGoAnnotationReference(
			Set<String> uniqueShortAccessions) throws IOException {
		Map<String, Set<String>> goa = new HashMap<String, Set<String>>();
		BufferedReader goaIn = null;
		for (String blastDatabaseName : getSettings().getBlastDatabases()) {
			if (getSettings().hasGeneOntologyAnnotation(blastDatabaseName)) {
				try {
					goaIn = new BufferedReader(new FileReader(getSettings()
							.getPathToGeneOntologyReference(blastDatabaseName)));
					Pattern p = getSettings().getGoReferenceRegex(blastDatabaseName);
					String line, shortAcc, goTerm = "";
					while ((line = goaIn.readLine()) != null) {
						Matcher m = p.matcher(line);
						if (m.find()) {
							shortAcc = m.group(SHORT_ACCESSION_GROUP_NAME);
							if (uniqueShortAccessions.contains(shortAcc)) {
								goTerm = m.group(GO_TERM_GROUP_NAME);
								addGoAnnotation(goa, shortAcc, goTerm);
							}
						}
					}
				} finally {
					goaIn.close();
				}
			}
		}
		return goa;
	}

	/**
	 * Adds the Gene Ontology term <code>goTerm</code> to the Set of the
	 * BlastResult's GO term annotations. In this, the BlastResult is identified
	 * by its short accession <code>brShortAccession</code>.
	 * 
	 * @param goa
	 * @param brShortAccession
	 * @param goTerm
	 */
	protected static void addGoAnnotation(Map<String, Set<String>> goa,
			String brShortAccession, String goTerm) {
		if (!goa.containsKey(brShortAccession)) {
			goa.put(brShortAccession, new HashSet<String>());
		}
		goa.get(brShortAccession).add(goTerm);
	}
}
