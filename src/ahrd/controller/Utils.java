package ahrd.controller;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Provides globally used utility-methods. E.g. for reading files or creating
 * Lists of Strings representing the lines of a File.
 * 
 * @author hallab, klee
 */
public class Utils {

	public static final Random random = new Random();

	/**
	 * Uses (double) Math.round(value * 100000) / 100000 to round to 5 digits
	 * after decimal point.
	 * 
	 * @param toBeFormatted
	 * @param nDecimalPlaces
	 * @return rounded Double
	 */
	public static Double roundToNDecimalPlaces(Double toBeFormatted,
			int nDecimalPlaces) {
		double decPlacesFact = Math.pow(10, nDecimalPlaces);
		return (double) Math.round(toBeFormatted * decPlacesFact)
				/ decPlacesFact;
	}

	public static String readFile(String path) throws IOException {
		FileInputStream stream = new FileInputStream(new File(path));
		try {
			FileChannel fc = stream.getChannel();
			MappedByteBuffer bb = fc.map(FileChannel.MapMode.READ_ONLY, 0,
					fc.size());
			/* Instead of using default, pass in a decoder. */
			return Charset.defaultCharset().decode(bb).toString();
		} finally {
			stream.close();
		}
	}

	/**
	 * Random: >= 0.1 and <= 1.0
	 * 
	 * @return Long
	 */
	public static Double randomMultipleOfOneTenth() {
		return randomMultipleOfTen() * 0.01;
	}

	/**
	 * Random: >= 10 and <= 100
	 * 
	 * @return Long
	 */
	public static Long randomMultipleOfTen() {
		Random rand = Utils.random;
		return new Long((rand.nextInt(10) + 1) * 10);
	}

	public static boolean randomTrueOrFalse() {
		Random rand = Utils.random;
		;
		return rand.nextBoolean();
	}

	public static boolean randomSaveSubtract(int from, int subtr) {
		return (from - subtr > 0 ? randomTrueOrFalse() : false);
	}

	public static boolean randomSaveSubtract(double from, double subtr) {
		return (from - subtr > 0 ? randomTrueOrFalse() : false);
	}

	/**
	 * Reads a Files content and returns a List of lines, discarding empty
	 * lines.
	 * 
	 * @param path2File
	 * @return List<String> of lines
	 * @throws IOException
	 */
	public static List<String> fromFile(String path2File) throws IOException {
		List<String> fromFile = new ArrayList<String>();
		String fileContent = readFile(path2File);
		for (String line : fileContent.split("\n")) {
			String lineContent = line.trim();
			if (lineContent != null && !lineContent.equals(""))
				fromFile.add(lineContent);
		}
		return fromFile;
	}

}
