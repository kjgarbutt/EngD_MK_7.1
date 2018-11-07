package sim;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;

import agents.Agent;

import com.linuxense.javadbf.*;

import java.io.*;

import sim.engine.SimState;
import sim.field.geo.GeomVectorField;
import sim.io.geo.ShapeFileImporter;
import sim.util.Bag;
import sim.util.geo.GeomPlanarGraph;
import sim.util.geo.GeomPlanarGraphEdge;
import sim.util.geo.MasonGeometry;
import sim.util.geo.AttributeValue;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.planargraph.Node;
import network.GeoNode;
import network.ListEdge;

import ec.util.MersenneTwisterFast;

/**
 *
 * "MK_7_" is iteration 7.1 of my EngD project model. It varies little from
 * MK_7, and just has updated details to include Gloucestershire date. It is
 * adapted from the MASON demo, "Gridlock", made by Sarah Wise, Mark Coletti,
 * and Andrew Crooks.
 * 
 * The model reads a number of GIS shapefiles and displays a road network, two
 * Environment Agency flood maps and a bespoke Open Source Vulnerability Index
 * (OSVI). The model reads in a .CSV and generates a predetermined number of
 * agents with set characteristics. The agents are placed on the road network
 * and are located at a Red Cross office. The model reads a separate .CSV and
 * assigns goal locations to each agent at random from a predetermined list. The
 * agents are assigned speeds at random. Once the model is started, the agents
 * move from A to B, then they change direction and head back to their start
 * position. The process repeats until the user quits.
 *
 * @author KJGarbutt
 *
 */
public class MK_7_1 extends SimState {

	//////////////////////////////////////////////////////////////////////////////
	///////////////////////////// MODEL PARAMETERS ///////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	private static final long serialVersionUID = -4554882816749973618L;
	public static double resolution = 5;// the granularity of the simulation
	// (fiddle around with this to merge nodes into one another)

	///////////////////////////// Containers /////////////////////////////////////
	public GeomVectorField baseLayer = new GeomVectorField();
	public GeomVectorField roads = new GeomVectorField();
	public GeomVectorField flood3 = new GeomVectorField();
	public GeomVectorField flood2 = new GeomVectorField();
	public GeomVectorField agentsLayer = new GeomVectorField();

	///////////////////////////// Network ////////////////////////////////////////
	public GeomPlanarGraph network = new GeomPlanarGraph();
	// Stores road network connections
	public GeomVectorField junctions = new GeomVectorField();
	// Stores nodes for road intersections
	HashMap<Integer, GeomPlanarGraphEdge> idsToEdges = new HashMap<Integer, GeomPlanarGraphEdge>();
	public HashMap<GeomPlanarGraphEdge, ArrayList<agents.Agent>> edgeTraffic = new HashMap<GeomPlanarGraphEdge, ArrayList<agents.Agent>>();
	public GeomVectorField mainagents = new GeomVectorField();

	///////////////////////////// Objects ////////////////////////////////////////
	public GeometryFactory fa = new GeometryFactory();

	long mySeed = 0;

	Envelope MBR = null;

	// Model ArrayLists for agents and OSVI Polygons
	ArrayList<agents.Agent> agentList = new ArrayList<agents.Agent>();
	ArrayList<Integer> assignedWards = new ArrayList<Integer>();
	ArrayList<Integer> visitedWards = new ArrayList<Integer>(); // TODO record visited LSOAs
	ArrayList<Polygon> polys = new ArrayList<Polygon>();
	ArrayList<String> csvData = new ArrayList<String>();

	// Here we force the agents to go to or from work at any time
	public boolean goToLSOA = true;

	public int activeCount;
	public int numActive = 0;

	public boolean getGoToLSOA() {
		return goToLSOA;
	}

	/*
	 * //Need to actually utilise this somewhere public static double
	 * temporalResolution_minutesPerTick = 1;//5; // minute per tick public double
	 * param_defaultSpeed = 200 * MK_7_1.temporalResolution_minutesPerTick; public
	 * double param_topSpeed = 1000 * MK_7_1.temporalResolution_minutesPerTick;
	 * public int param_reportTimeCommitment = (int)(60 /
	 * temporalResolution_minutesPerTick); public int
	 * param_responseCarTimeCommitment = (int)(60 /
	 * temporalResolution_minutesPerTick);
	 */

	//////////////////////////////////////////////////////////////////////////////
	/////////////////////////// BEGIN FUNCTIONS //////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	/**
	 * //////////////////////// Model Constructor ////////////////////////////////
	 * Model Constructor
	 */
	public MK_7_1(long seed) {
		super(seed);
		random = new MersenneTwisterFast(12345);
	}

	/**
	 * //////////////////////// OSVI Polygon Setup ///////////////////////////////
	 * Polygon Setup
	 */
	void setup() {
		// copy over the geometries into a list of Polygons
		Bag ps = baseLayer.getGeometries();
		polys.addAll(ps);
	}

	/**
	 * //////////////////////// Model Initialisation /////////////////////////////
	 * Model Initialisation
	 */
	@Override
	public void start() {
		super.start();

		System.out.println();
		System.out.println("////////////////\nINPUTTING STUFFS\n////////////////");
		System.out.println();

		System.out.println("Reading shapefiles...");

		//////////////////////////////////////////////////////////////////////////
		/////////////////////////// READ IN DATA /////////////////////////////////
		//////////////////////////////////////////////////////////////////////////

		try {
			File wardsFile = new File("data/GloucestershireFinal_LSOA1.shp");
			ShapeFileImporter.read(wardsFile.toURI().toURL(), baseLayer, Polygon.class);
			System.out.println("	OSVI shapefile: " + wardsFile);

			// read in the roads shapefile to create the transit network
			File roadsFile = new File("data/GL_ITN_MultipartToSinglepart.shp");
			ShapeFileImporter.read(roadsFile.toURI().toURL(), roads);
			System.out.println("	Roads shapefile: " + roadsFile);

			/*
			 * // read in the LSOA shapefile to create the backgrounds URL areasFile =
			 * MK_4.class.getResource ("/data/Final_LSOA.shp"); Bag desiredAttributes = new
			 * Bag(); desiredAttributes.add("RC_RankCol"); try {
			 * ShapeFileImporter.read(areasFile, world, desiredAttributes); } catch
			 * (FileNotFoundException ex){ }
			 */

			// read in the FZ3 file
			File flood3File = new File("data/Gloucestershire_FZ_3.shp");
			ShapeFileImporter.read(flood3File.toURI().toURL(), flood3);
			System.out.println("	FZ3 shapefile: " + flood3File);

			// read in the FZ2 file
			File flood2File = new File("data/Gloucestershire_FZ_2.shp");
			ShapeFileImporter.read(flood2File.toURI().toURL(), flood2);
			System.out.println("	FZ2 shapefile: " + flood2File);

			createNetwork();

			System.out.println("Setting up OSVI Portrayals...");
			System.out.println();

			setup();

			//////////////////////////////////////////////////////////////////////
			/////////////////////////// CLEANUP //////////////////////////////////
			//////////////////////////////////////////////////////////////////////

			// clear any existing agents from previous runs
			agentsLayer.clear();

			// standardize the MBRs so that the visualization lines up
			// and everyone knows what the standard MBR is
			MBR = baseLayer.getMBR();
			baseLayer.setMBR(MBR);
			roads.setMBR(MBR);
			flood3.setMBR(MBR);
			flood2.setMBR(MBR);

			//////////////////////////////////////////////////////////////////////
			/////////////////////////// AGENTS ///////////////////////////////////
			//////////////////////////////////////////////////////////////////////

			// initialize agents using the following source .CSV files
			setupAgentsFromFile("data/GloucestershireITNAGENT.csv");
			agentsLayer.setMBR(MBR);

			System.out.println();
			System.out.println("Starting simulation...");

			// Ensure that the spatial index is updated after all the agents move
			schedule.scheduleRepeating(agentsLayer.scheduleSpatialIndexUpdater(), Integer.MAX_VALUE, 1.0);

			/**
			 * Steppable that flips Agent paths once everyone reaches their destinations
			 * 
			 * Steppable flipper = new Steppable() {
			 * 
			 * @Override public void step(SimState state) {
			 * 
			 *           MK_7 gstate = (MK_7) state;
			 * 
			 *           // pass to check if anyone has not yet reached work //for
			 *           (MainAgent a : gstate.agentList) //{ ///if (!a.reachedDestination)
			 *           //{ // return; // someone is still moving: let him do so // } //}
			 *           // send everyone back in the opposite direction now //boolean
			 *           toLSOA = gstate.goToLSOA; // gstate.goToLSOA = !toLSOA;
			 * 
			 *           // otherwise everyone has reached their latest destination: // turn
			 *           them back for (MainAgent a : gstate.agentList) if
			 *           (a.reachedDestination) { a.flipPath(); } } };
			 *           schedule.scheduleRepeating(flipper, 10);
			 */

		} catch (FileNotFoundException e) {
			System.out.println("Error: missing required data file");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * //////////////////////// Model Finish & Cleanup ///////////////////////////
	 * Finish the simulation and clean up
	 */
	public void finish() {
		super.finish();

		System.out.println();
		System.out.println("Simulation ended by user.");

		System.out.println();
		System.out.println("///////////////////////\nOUTPUTTING FINAL STUFFS\n///////////////////////");
		System.out.println();

		for (Agent a : agentList) {
			for (String s : a.recordOfTrips)
				System.out.println(s);
		}
	}

	/**
	 * //////////////////////// Create Road Network //////////////////////////////
	 * Create the road network the agents will traverse
	 */
	private void createNetwork() {
		System.out.println("Creating road network..." + roads);
		System.out.println();
		network.createFromGeomField(roads);

		for (Object o : network.getEdges()) {
			GeomPlanarGraphEdge e = (GeomPlanarGraphEdge) o;

			idsToEdges.put(e.getIntegerAttribute("ROAD_ID_1").intValue(), e);

			e.setData(new ArrayList<agents.Agent>());
		}

		addIntersectionNodes(network.nodeIterator(), junctions);
	}

	/**
	 * ///////////////////////// Setup agentGoals /////////////////////////////////
	 * Read in the agent goals CSV
	 * 
	 * @param agentfilename
	 * @return
	 *
	 */
	public ArrayList<String> agentGoals(String agentfilename) throws IOException {
		String csvGoal = null;
		BufferedReader agentGoalsBuffer = null;

		String agentFilePath = MK_7_1.class.getResource(agentfilename).getPath();
		FileInputStream agentfstream = new FileInputStream(agentFilePath);
		System.out.println("Reading Agent's Goals file: " + agentFilePath);

		try {
			agentGoalsBuffer = new BufferedReader(new InputStreamReader(agentfstream));
			agentGoalsBuffer.readLine();
			while ((csvGoal = agentGoalsBuffer.readLine()) != null) {
				String[] splitted = csvGoal.split(",");

				ArrayList<String> agentGoalsResult = new ArrayList<String>(splitted.length);
				for (String data : splitted)
					agentGoalsResult.add(data);
				csvData.addAll(agentGoalsResult);
			}
			System.out.println();
			System.out.println("Full csvData Array: " + csvData);

		} finally {
			if (agentGoalsBuffer != null)
				agentGoalsBuffer.close();
		}
		return csvData;
	}

	int getLargestUnassignedWard() {
		Bag lsoaGeoms = baseLayer.getGeometries();

		int highestOSVI = -1;
		MasonGeometry myCopy = null;

		for (Object o : lsoaGeoms) {
			MasonGeometry masonGeometry = (MasonGeometry) o;
			int id = masonGeometry.getIntegerAttribute("ID");
			if (assignedWards.contains(id))
				continue;

			int tempOSVI = masonGeometry.getIntegerAttribute("L_GL_OSVI_");
			// temp = the attribute in the "L_GL_OSVI_" column (int for each LSOA OSVI)
			if (tempOSVI > highestOSVI) { // if temp is higher than highest
				highestOSVI = tempOSVI; // update highest to temp
				myCopy = masonGeometry; // update myCopy to the
			}

			if (myCopy == null) {
				System.out.println("ALERT: LSOA Baselayer is null!");
			}
		}
		
		///////////////////////////////////////////////////////////
		////// TODO HOW TO STOP myCopy ENDING UP AT NULL??? ///////
		///////////////////////////////////////////////////////////

		int id = myCopy.getIntegerAttribute("ID");
		// id = the attribute in the "ID" column (int for each LSOA)
		assignedWards.add(id); // add ID to the "assignedWards" ArrayList
		System.out.println();
		System.out.println("Highest OSVI Raiting: " + myCopy.getIntegerAttribute("L_GL_OSVI_"));
		return myCopy.getIntegerAttribute("ROAD_ID_1"); // return Road_ID for the chosen LSOA to visit
	}

	/**
	 * //////////////////////// Setup Agents ////////////////////////////////// Read
	 * in the population files and create appropriate populations
	 * 
	 * @param agentsFilename
	 */
	public void setupAgentsFromFile(String agentsFilename) {
		try {
			// String filePath = MK_7_1.class.getResource(filename).getPath();
			FileInputStream fstream = new FileInputStream(agentsFilename);
			System.out.println("Reading in Agents from: " + agentsFilename);

			BufferedReader d = new BufferedReader(new InputStreamReader(fstream));
			String s;

			// get rid of the header
			d.readLine();
			// read in all data
			while ((s = d.readLine()) != null) {
				String[] bits = s.split(",");

				int pop = Integer.parseInt(bits[0]);

				// moveRate = (int)(Math.random()*70) + 1;
				// System. out.println("MoveRate = " + moveRate );
				// int mainAgentSpeed = MainAgent.MoveRate;
				// System.out.println("Main Agent speed = " +mainAgentSpeed);

				String id = bits[3];
				String homeTract = bits[1];
				String ROAD_ID = bits[1];
				String agentName = bits[2];

				Random randomiser = new Random();
				int random = getLargestUnassignedWard();
				// String random = csvData.get(new Random().nextInt(csvData.size()));
				// String goalTract = random;
				System.out.println("Nearest road segment of of chosen LSOA Centroid: " + random);// goalTract);
				System.out.println("........");
				System.out.println("Assigning Agent to chosen LSOA... ");
				System.out.println("Agent's name: " + agentName);// goalTract);
				System.out.println("Agent's HQ road segment: " + homeTract);

				GeomPlanarGraphEdge startingEdge = idsToEdges.get((int) Double.parseDouble(ROAD_ID));
				GeomPlanarGraphEdge goalEdge = idsToEdges.get(random);// (int) Double.parseDouble(goalTract));
				// reads the .CSV column
				// goals[ random.nextInt(goals.length)]);
				// uses the hardcoded 'goals' from above

				for (int i = 0; i < pop; i++) {
					agents.Agent a = new agents.Agent(agentName, null, agentName, startingEdge, goalEdge);

					boolean successfulStart = a.start(this);
					// System.out.println("Starting...");

					if (!successfulStart) {
						System.out.println("ERROR: Main agents *NOT* added properly!");
						continue; // DON'T ADD IT if it's bad
					} else {
						// System.out.println("Agent added successfully!");
					}

					// MasonGeometry newGeometry = new MasonGeometry(a.getGeometry());
					MasonGeometry newGeometry = a.getGeometry();
					newGeometry.isMovable = true;
					agentsLayer.addGeometry(newGeometry);
					agentList.add(a);
					schedule.scheduleOnce(a);
				}
			}

			d.close();
			System.out.println();
			System.out.println("All Agents added successfully!");
		} catch (Exception e) {
			System.out.println("ERROR: issue with population file: ");
			e.printStackTrace();
		}
	}

	/**
	 * //////////////////////// Network Intersections ////////////////////////////
	 * adds nodes corresponding to road intersections to GeomVectorField
	 *
	 * @param nodeIterator
	 *            Points to first node
	 * @param intersections
	 *            GeomVectorField containing intersection geometry
	 *
	 *            Nodes will belong to a planar graph populated from LineString
	 *            network.
	 */
	private void addIntersectionNodes(Iterator<?> nodeIterator, GeomVectorField intersections) {
		System.out.println("Adding Intersection Nodes...");
		System.out.println();
		GeometryFactory fact = new GeometryFactory();
		Coordinate coord = null;
		Point point = null;
		@SuppressWarnings("unused")
		int counter = 0;

		while (nodeIterator.hasNext()) {
			Node node = (Node) nodeIterator.next();
			coord = node.getCoordinate();
			point = fact.createPoint(coord);

			junctions.addGeometry(new MasonGeometry(point));
			counter++;
		}
	}

	/**
	 * Return the GeoNode in the road network which is closest to the given
	 * coordinate
	 * 
	 * @param c
	 * @return
	 */
	public GeoNode getClosestGeoNode(Coordinate c) {

		// find the set of all nodes within *resolution* of the given point
		Bag objects = roads.getObjectsWithinDistance(fa.createPoint(c), resolution);
		if (objects == null || roads.getGeometries().size() <= 0)
			return null; // problem with the network layer

		// among these options, pick the best
		double bestDist = resolution; // MUST be within resolution to count
		GeoNode best = null;
		for (Object o : objects) {
			double dist = ((GeoNode) o).geometry.getCoordinate().distance(c);
			if (dist < bestDist) {
				bestDist = dist;
				best = ((GeoNode) o);
			}
		}

		// if there is a best option, return that!
		if (best != null && bestDist == 0)
			return best;

		// otherwise, closest GeoNode is associated with the closest Edge, so look for
		// that!

		ListEdge edge = getClosestEdge(c);

		// find that edge
		if (edge == null) {
			edge = getClosestEdge(c, resolution * 10);
			if (edge == null)
				return null;
		}

		// of that edge's endpoints, find the closer of the two and return it
		GeoNode n1 = (GeoNode) edge.getFrom();
		GeoNode n2 = (GeoNode) edge.getTo();

		if (n1.geometry.getCoordinate().distance(c) <= n2.geometry.getCoordinate().distance(c))
			return n1;
		else
			return n2;
	}

	/**
	 * Return the ListEdge in the road network which is closest to the given
	 * coordinate
	 * 
	 * @param c
	 * @return
	 */
	public ListEdge getClosestEdge(Coordinate c) {

		// find the set of all edges within *resolution* of the given point
		Bag objects = roads.getObjectsWithinDistance(fa.createPoint(c), resolution);
		if (objects == null || roads.getGeometries().size() <= 0)
			return null; // problem with the network edge layer

		Point point = fa.createPoint(c);

		// find the closest edge among the set of edges
		double bestDist = resolution;
		ListEdge bestEdge = null;
		for (Object o : objects) {
			double dist = ((MasonGeometry) o).getGeometry().distance(point);
			if (dist < bestDist) {
				bestDist = dist;
				bestEdge = (ListEdge) ((AttributeValue) ((MasonGeometry) o).getAttribute("ListEdge")).getValue();
			}
		}

		// if it exists, return it
		if (bestEdge != null)
			return bestEdge;

		// otherwise return failure
		else
			return null;
	}

	/**
	 * Return the ListEdge in the road network which is closest to the given
	 * coordinate, within the given resolution
	 * 
	 * @param c
	 * @param resolution
	 * @return
	 */
	public ListEdge getClosestEdge(Coordinate c, double resolution) {

		// find the set of all edges within *resolution* of the given point
		Bag objects = roads.getObjectsWithinDistance(fa.createPoint(c), resolution);
		if (objects == null || roads.getGeometries().size() <= 0)
			return null; // problem with the network edge layer

		Point point = fa.createPoint(c);

		// find the closest edge among the set of edges
		double bestDist = resolution;
		ListEdge bestEdge = null;
		for (Object o : objects) {
			double dist = ((MasonGeometry) o).getGeometry().distance(point);
			if (dist < bestDist) {
				bestDist = dist;
				bestEdge = (ListEdge) ((AttributeValue) ((MasonGeometry) o).getAttribute("ListEdge")).getValue();
			}
		}

		// if it exists, return it
		if (bestEdge != null)
			return bestEdge;

		// otherwise return failure
		else
			return null;
	}

	/** set the seed of the random number generator */
	void seedRandom(long number) {
		random = new MersenneTwisterFast(number);
		mySeed = number;
	}

	/**
	 * //////////////////////// Main Function ////////////////////////////////////
	 * Main function allows simulation to be run in stand-alone, non-GUI mode
	 */
	public static void main(String[] args) {
		doLoop(MK_7_1.class, args);
		System.exit(0);
	}
}
