package ahrd.model;

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

public class ReferenceGoAnnotation {

	private static final String SHORT_ACCESSION_GROUP_NAME = "shortAccession";
	private static final String GO_TERM_GROUP_NAME = "goTerm";
	private static final String EVIDENCE_CODE_GROUP_NAME = "evidenceCode";
	private String goTerm;
	private String evidenceCode;

	public ReferenceGoAnnotation(String term, String code) {
		super();
		this.goTerm = term;
		this.evidenceCode = code;
	}
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
	public static Map<String, Set<ReferenceGoAnnotation>> parseGoAnnotationReference(
			Set<String> uniqueShortAccessions) throws IOException {
		Map<String, Set<ReferenceGoAnnotation>> goa = new HashMap<String, Set<ReferenceGoAnnotation>>();
		BufferedReader goaIn = null;
		for (String blastDatabaseName : getSettings().getBlastDatabases()) {
			if (getSettings().hasGeneOntologyAnnotation(blastDatabaseName)) {
				try {
					goaIn = new BufferedReader(new FileReader(getSettings()
							.getPathToGeneOntologyReference(blastDatabaseName)));
					Pattern p = getSettings().getGoReferenceRegex(blastDatabaseName);
					String line, shortAcc, term, code = "";
					while ((line = goaIn.readLine()) != null) {
						Matcher m = p.matcher(line);
						if (m.find()) {
							shortAcc = m.group(SHORT_ACCESSION_GROUP_NAME);
							if (uniqueShortAccessions.contains(shortAcc)) {
								term = m.group(GO_TERM_GROUP_NAME);
								code = m.group(EVIDENCE_CODE_GROUP_NAME);
								if (!goa.containsKey(shortAcc)) {
									goa.put(shortAcc, new HashSet<ReferenceGoAnnotation>());
								}
								goa.get(shortAcc).add(new ReferenceGoAnnotation(term, code));
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

	public String getGoTerm() {
		return goTerm;
	}

	public void setGoTerm(String goTerm) {
		this.goTerm = goTerm;
	}

	public String getEvidenceCode() {
		return evidenceCode;
	}

	public void setEvidenceCode(String evidenceCode) {
		this.evidenceCode = evidenceCode;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj == this) {
			return true;
		}
		if((obj != null) && (obj.getClass() == this.getClass())) {
			if (this.goTerm.equals(((ReferenceGoAnnotation)obj).getGoTerm())) {
				return true;
			}
		}
		return false;
	}
	
	@Override
    public int hashCode() {
		return this.goTerm.hashCode();
	}
}
