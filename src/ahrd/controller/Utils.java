package ahrd.controller;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.Vector;

import nu.xom.Attribute;
import nu.xom.Element;
import nu.xom.Nodes;
import nu.xom.XPathContext;

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

	public static List<Double> zeroList(int size) {
		List<Double> l = new Vector<Double>(size);
		for (int i = 0; i < size; i++) {
			l.add(0.0);
		}
		return l;
	}

	/**
	 * Returns the value XML attribute 'attributeName' of 'element'. Returns
	 * EMPTY String in case the attribute can't be found.
	 * 
	 * @param element
	 * @param attributeName
	 * @return String
	 */
	public static String getXmlAttributeValue(Element element,
			String attributeName) {
		String attrVal = "";
		Attribute attr = element.getAttribute(attributeName);
		if (attr != null)
			attrVal = attr.getValue();
		return attrVal;
	}

	/**
	 * Returns the content of the <em>first</em> XML child of the argument
	 * element.
	 * 
	 * @param element
	 * @param xpathQuery
	 * @return String
	 */
	public static String retrieveContentOfFirstXmlChildElement(Element element,
			String xpathQuery) {
		String res = null;
		Nodes resultNodes = element.query(xpathQuery);
		if (resultNodes.size() > 0) {
			res = resultNodes.get(0).getValue();
		}
		return res;
	}

	/**
	 * Returns the content of the <em>first</em> XML child of the argument
	 * element. Uses provided argument name space.
	 * 
	 * @param element
	 * @param xpathQuery
	 * @param context
	 * @return String
	 */
	public static String retrieveContentOfFirstXmlChildElement(Element element,
			String xpathQuery, XPathContext context) {
		String res = null;
		Nodes resultNodes = element.query(xpathQuery, context);
		if (resultNodes.size() > 0) {
			res = resultNodes.get(0).getValue();
		}
		return res;
	}

	/**
	 * Returns all values of attributes matching argument name found in the
	 * children of argument XML element.
	 * 
	 * @param element
	 * @param xpathQuery
	 * @param attributeName
	 * @return Set<String>
	 */
	public static Set<String> retrieveAttribteValuesOfXmlChildrenElements(
			Element element, String xpathQuery, String attributeName) {
		Set<String> res = null;
		Nodes resultNodes = element.query(xpathQuery);
		if (resultNodes.size() > 0) {
			res = new HashSet<String>();
			for (int i = 0; i < resultNodes.size(); i++) {
				res.add(getXmlAttributeValue((Element) resultNodes.get(i),
						attributeName));
			}
		}
		return res;
	}

	/**
	 * Returns all values of attributes matching argument name found in the
	 * children of argument XML element. Uses provided argument name space.
	 * 
	 * @param element
	 * @param xpathQuery
	 * @param attributeName
	 * @param context
	 * @return Set<String>
	 */
	public static Set<String> retrieveAttribteValuesOfXmlChildrenElements(
			Element element, String xpathQuery, String attributeName,
			XPathContext context) {
		Set<String> res = null;
		Nodes resultNodes = element.query(xpathQuery, context);
		if (resultNodes.size() > 0) {
			res = new HashSet<String>();
			for (int i = 0; i < resultNodes.size(); i++) {
				res.add(getXmlAttributeValue((Element) resultNodes.get(i),
						attributeName));
			}
		}
		return res;
	}
}
