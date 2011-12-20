package barna.astalavista;

import barna.model.splicegraph.SplicingGraph;

import java.lang.reflect.Method;
import java.util.HashMap;


public class AStalavista {

	public static final String VERSION_ID= "2.2";
	public static final Class[] mainSig= new Class[] {String[].class};
	public static HashMap<String, Class> commandMap=
		new HashMap<String, Class>(1,1f);
	public final static String commandValidStr= "asta lavista sort extractSJ extractAttributes";
	// TODO took out GFFSorter.class to make it work with current code
	public final static Class[] commandClass= new Class[] {SplicingGraph.class, null, SJextractor.class, AttributeExtractor.class};	// null= LaVista.class
	public static String commandDescription= "The command to be executed. ";
	static {
		
		String[] commandValid= commandValidStr.split("\\s");
		for (int i = 0; i < commandValid.length; i++) {
			commandMap.put(commandValid[i], commandClass[i]);
			commandDescription+= commandValid[i];
			if (i+1< commandValid.length)
				commandDescription+= ", ";
		}
		
	}
	
	public static String[] addHelp(String[] args) {
		 String[] newArgs= new String[args.length+1];
		 System.arraycopy(args, 0, newArgs, 0, args.length);
		 newArgs[newArgs.length-1]= "--help";
		 return newArgs;
	}
	
	public static void main(String[] args) {
		
		System.err.println("Welcome to the AStalavista toolkit! (Version "+VERSION_ID+")");
		System.err.print(AStalavista.class.getName()+" ");
		for (int i = 0; i < args.length; i++) 
			System.err.print(args[i]+" ");
		System.err.println();
		Class c= null;
		if (args!= null&& args.length> 1)
			c= commandMap.get(args[1]);
//		for (int i = 0; i < args.length; i++) {
//			if (args[i].equals("-c")|| args[i].equals("--command")&& (i+1)<  args.length) {
//				System.err.println("command "+args[i+1]);
//				c= commandMap.get(args[i+1]);
//				++i;
//				break;
//			}
//		} 
		
		if (c== null) {
			System.err.println("Usage: -c --command value: asta");
			System.exit(0);
		}
		
		try {			
			Method m= getMainMethod(c);
			// String[] rest= opt.getUnparsedArgs(); // geht leider net
			String[] args1= new String[args.length- 2];
			System.arraycopy(args, 2, args1, 0, args1.length);
			m.invoke(null, new Object[] {args1});
		} catch (Exception ex) {
			ex.printStackTrace();
			System.err.println("astalavista.");
			System.exit(-1);
		}
		 
	}

	static Method getMainMethod(Class c) throws SecurityException, NoSuchMethodException {
		return c.getDeclaredMethod("main", mainSig);
	}
}
