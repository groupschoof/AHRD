package ahrd.controller;

import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.util.List;

import ahrd.model.Protein;
import ahrd.model.ReferenceDescription;

public class FastaIndexer {

	public static List<ReferenceDescription> references;
	
	public static void main(String[] args) {
		List<String> fastaEntries = Protein.splitFasta(args[0]);
		for (String fastaEntry : fastaEntries) {
			if (fastaEntry != null && !fastaEntry.trim().equals("")) {
				references.add(ReferenceDescription
						.constructFromFastaEntry(fastaEntry.trim()));
			}
		}
		try {
			ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(args[1]));
			out.writeObject(references);
			out.close();
		} catch (Exception e) {
			e.printStackTrace(System.err);
		}
	}

}
