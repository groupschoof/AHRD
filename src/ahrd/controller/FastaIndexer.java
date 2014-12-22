package ahrd.controller;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Properties;

public class FastaIndexer {

	public static long timestamp;

	public static void main(String[] args) throws Exception {
		timestamp = (new Date()).getTime();
		System.out.println("Start");
		BufferedReader fastaIn = new BufferedReader(new FileReader(args[0]));
		String str;
		PrintWriter pw = new PrintWriter(new FileWriter(args[1]));
		while ((str = fastaIn.readLine()) != null) {
			if (str.startsWith(">")) {
				pw.println(str.trim().replaceFirst(">", "")
						.replaceFirst(" ", "="));
			}
		}
		fastaIn.close();
		pw.close();
		System.out.println("Serialized in " + takeTime());

		Properties prop = new Properties();
		InputStream stream = new FileInputStream(args[1]);
		prop.load(stream);
		System.out.println("Deserialized in " + takeTime());

		// System.out.println(references);
	}

	static protected long takeTime() {
		// Measure time:
		long now = (new Date()).getTime();
		long measuredSeconds = (now - timestamp) / 1000;
		timestamp = now;
		return measuredSeconds;
	}

}
