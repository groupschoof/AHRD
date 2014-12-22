package ahrd.controller;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ahrd.model.ReferenceDescription;

public class FastaIndexer {

	public static List<Map<String, String>> references = new ArrayList<Map<String, String>>();
	public static long timestamp;

	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws Exception {
		timestamp = (new Date()).getTime();
		System.out.println("Start");
		BufferedReader fastaIn = new BufferedReader(new FileReader(args[0]));
		String str;
		while ((str = fastaIn.readLine()) != null) {
			if (str.startsWith(">")) {
				ReferenceDescription rd = ReferenceDescription
						.constructFromFastaEntry(str.trim());
				Map<String, String> m = new HashMap<String, String>();
				m.put(rd.getAccession(), rd.getDescription());
				references.add(m);
			}
		}
		fastaIn.close();
		ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(
				args[1]));
		out.writeObject(references);
		out.close();
		System.out.println("Serialized in " + takeTime());

		references = new ArrayList<Map<String, String>>();
		ObjectInputStream in = new ObjectInputStream(new FileInputStream(
				args[1]));
		references = (List<Map<String, String>>) in.readObject();
		in.close();
		System.out.println("Deserialized in " + takeTime());
		
		System.out.println(references);
	}
	
	static protected long takeTime() {
		// Measure time:
		long now = (new Date()).getTime();
		long measuredSeconds = (now - timestamp) / 1000;
		timestamp = now;
		return measuredSeconds;
	}

}
