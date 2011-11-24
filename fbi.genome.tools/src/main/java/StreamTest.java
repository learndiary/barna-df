
public class StreamTest {

	public static void main(String[] args) {
		for (int i = 0; i < 10; i++) {
			System.err.print("working [");
			for (int j = 0; j < i; j++) 
				System.err.print('-');
			System.err.print('\\');
			for (int j = i+1; j < 10; j++) 
				System.err.print(' ');
			System.err.print("]\r");
			System.err.flush();
			
			try {
				Thread.currentThread().sleep(100);
			} catch (InterruptedException e) {
				; // :)
			}

			System.err.print("working [");
			for (int j = 0; j < i; j++)
				System.err.print('-');
			System.err.print('|');
			for (int j = i+1; j < 10; j++) 
				System.err.print(' ');
			System.err.print("]\r");
			System.err.flush();
			
			try {
				Thread.currentThread().sleep(100);
			} catch (InterruptedException e) {
				; // :)
			}

			System.err.print("working [");
			for (int j = 0; j < i; j++)
				System.err.print('-');
			System.err.print('/');
			for (int j = i+1; j < 10; j++) 
				System.err.print(' ');
			System.err.print("]\r");
			System.err.flush();
			
			try {
				Thread.currentThread().sleep(100);
			} catch (InterruptedException e) {
				; // :)
			}

			System.err.print("working [");
			for (int j = 0; j <= i; j++) 
				System.err.print('-');
			for (int j = i+1; j < 10; j++) 
				System.err.print(' ');
			System.err.print("]\r");
			System.err.flush();
			
			try {
				Thread.currentThread().sleep(100);
			} catch (InterruptedException e) {
				; // :)
			}
			
		}
	}
}
