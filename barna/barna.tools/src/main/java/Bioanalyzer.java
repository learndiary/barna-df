import java.io.ByteArrayInputStream;
import java.io.EOFException;
import java.io.IOException;

import com.mindprod.ledatastream.LEDataInputStream;

public class Bioanalyzer {

	/**
	 * Expert Software XML Schema Date printed: 3/18/11 Agilent Technologies
	 * Page 15 of 54 Author: Volker von Einem
	 */
	static void originalSnippet() {
		
		String localName= null, content= null;
		
		if (localName.equals("RawSignal")) {
			
			ByteArrayInputStream bAIS = new ByteArrayInputStream(
					content.getBytes());
			
			Base64.InputStream b64IStream = new Base64.InputStream(bAIS,
					Base64.DECODE);
			
			LEDataInputStream lEDIStream = new LEDataInputStream(b64IStream);
			
			double yValue = 0;
			boolean endOfStream = false;
			while (!endOfStream) {
				try {
					// read values bit by bit
					yValue = (double) lEDIStream.readFloat();
					// use values
				} catch (EOFException e) {
					endOfStream = true;
				} catch (IOException e) {
					; // was uncaught in snippet
				}
			}
		}
	}
}
