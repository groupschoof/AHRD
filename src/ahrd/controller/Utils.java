package ahrd.controller;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
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
	 * Random: >= 1 and <= 10
	 * 
	 * @return Long
	 */
	public static Double randomMultipleOfOne() {
		return randomMultipleOfTen() * 0.1;
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
	
	/**
	 * From: http://stackoverflow.com/questions/15359702/get-location-of-jar-file/20953376#20953376
	 * 
	 * Compute the absolute file path to the jar file. The framework is based on
	 * http://stackoverflow.com/a/12733172/1614775 But that gets it right for
	 * only one of the four cases.
	 * 
	 * @param aclass
	 *            A class residing in the required jar.
	 * 
	 * @return A File object for the directory in which the jar file resides.
	 *         During testing with NetBeans, the result is ./build/classes/,
	 *         which is the directory containing what will be in the jar.
	 */
	@SuppressWarnings("rawtypes")
	public static String getAHRDdir(Class aclass) {
		URL url;
		String extURL; // url.toExternalForm();

		// get an url
		try {
			url = aclass.getProtectionDomain().getCodeSource().getLocation();
			// url is in one of two forms
			// ./build/classes/ NetBeans test
			// jardir/JarName.jar from a jar
		} catch (SecurityException ex) {
			url = aclass.getResource(aclass.getSimpleName() + ".class");
			// url is in one of two forms, both ending
			// "/com/physpics/tools/ui/PropNode.class"
			// file:/U:/Fred/java/Tools/UI/build/classes
			// jar:file:/U:/Fred/java/Tools/UI/dist/UI.jar!
		}

		// convert to external form
		extURL = url.toExternalForm();

		// prune for various cases
		if (extURL.endsWith(".jar")) // from getCodeSource
			extURL = extURL.substring(0, extURL.lastIndexOf("dist/"));
		else { // from getResource
			String suffix = (aclass.getName()).replace(".", "/") + ".class";
			extURL = replaceLast(extURL, suffix , "");
			//extURL = extURL.replace(suffix, "");
			suffix = "classes/";
			extURL = replaceLast(extURL, suffix , "");
			if (extURL.startsWith("jar:") && extURL.endsWith(".jar!"))
				extURL = extURL.substring(4, extURL.lastIndexOf("dist/"));
		}

		// convert back to url
		try {
			url = new URL(extURL);
		} catch (MalformedURLException mux) {
			// leave url unchanged; probably does not happen
		}

		return url.getPath();
	}
	
	public static String replaceLast(String text, String regex, String replacement) {
        return text.replaceFirst("(?s)(.*)" + regex, "$1" + replacement);
    }
	
	/**
	 * https://rosettacode.org/wiki/Evaluate_binomial_coefficients#Java
	 * (GNU Free Documentation License 1.2) 
	 */
	public static long binomial(int n, int k) {
        if (k>n-k)
            k=n-k;
 
        long b=1;
        for (int i=1, m=n; i<=k; i++, m--)
            b=b*m/i;
        return b;
    }
}
