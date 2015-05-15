package ahrd.view;

import static ahrd.controller.Settings.getSettings;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import ahrd.model.GOterm;
import ahrd.model.Protein;

public class ExtendedGOAnnotationTableWriter {

	public static final String HEADER = "Protein\tGO-Accession\tGO-Ontology\tGO-Name\n";

	private Collection<Protein> proteins;
	private Map<String, GOterm> goDb;

	public ExtendedGOAnnotationTableWriter(Collection<Protein> proteins,
			Map<String, GOterm> goDb) {
		super();
		setProteins(proteins);
		setGoDb(goDb);
	}

	public void writeOutput() throws IOException {
		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter(getSettings()
					.getExtendedGoResultTablePath()));
			bw.write(HEADER);
			for (Protein p : proteins) {
				for (String goAcc : p.getGoResults()) {
					GOterm g = getGoDb().get(goAcc);
					bw.write(p.getAccession() + "\t" + g.getAccession() + "\t"
							+ g.getOntology() + "\t" + g.getName() + "\n");
				}
			}
		} finally {
			bw.close();
		}
	}

	public Collection<Protein> getProteins() {
		return proteins;
	}

	public void setProteins(Collection<Protein> proteins) {
		this.proteins = proteins;
	}

	public Map<String, GOterm> getGoDb() {
		return goDb;
	}

	public void setGoDb(Map<String, GOterm> goDb) {
		this.goDb = goDb;
	}
}