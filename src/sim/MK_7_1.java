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
import sim.field.network.Edge;
import sim.field.network.Network;
import sim.io.geo.ShapeFileImporter;
import sim.util.Bag;
import sim.util.geo.GeomPlanarGraph;
import sim.util.geo.GeomPlanarGraphEdge;
import sim.util.geo.MasonGeometry;
import utilities.NetworkUtilities;
import sim.util.geo.AttributeValue;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.planargraph.Node;
import network.GeoNode;
import network.ListEdge;

import ec.util.MersenneTwisterFast;
import 	network.AStar;

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
	public int grid_width = 700;
	public int grid_height = 700;
	public static double resolution = 5;// the granularity of the simulation
	// (fiddle around with this to merge nodes into one another)

	///////////////////////////// Containers /////////////////////////////////////
	public GeomVectorField baseLayer = new GeomVectorField();
	public GeomVectorField roadLayer = new GeomVectorField();
	public GeomVectorField networkLayer = new GeomVectorField();
	public GeomVectorField networkEdgeLayer = new GeomVectorField();	
	public GeomVectorField majorRoadNodesLayer = new GeomVectorField();
	public GeomVectorField flood3 = new GeomVectorField();
	public GeomVectorField flood2 = new GeomVectorField();
	public GeomVectorField agentsLayer = new GeomVectorField();

	///////////////////////////// Network ////////////////////////////////////////
	//public GeomPlanarGraph network = new GeomPlanarGraph();
	// Stores road network connections
	public GeomVectorField junctions = new GeomVectorField();
	// Stores nodes for road intersections
	HashMap<Integer, GeomPlanarGraphEdge> idsToEdges = new HashMap<Integer, GeomPlanarGraphEdge>();
	public HashMap<GeomPlanarGraphEdge, ArrayList<agents.Agent>> edgeTraffic = new HashMap<GeomPlanarGraphEdge, ArrayList<agents.Agent>>();
	public GeomVectorField mainagents = new GeomVectorField();

	///////////////////////////// Objects ////////////////////////////////////////
	
	public Bag roadNodes = new Bag();
	public Network roadNetwork = new Network(false);
	HashMap <MasonGeometry, ArrayList <GeoNode>> localNodes;
	public Bag terminus_points = new Bag();
	
	public GeometryFactory fa = new GeometryFactory();

	long mySeed = 0;

	Envelope MBR = null;

	// Here we force the agents to go to or from work at any time
	public boolean goToLSOA = true;
	
	boolean verbose = false;

	public int activeCount;
	public int numActive = 0;
	
	// Model ArrayLists for agents and OSVI Polygons
	ArrayList<agents.Agent> agentList = new ArrayList<agents.Agent>();
	ArrayList<Integer> assignedWards = new ArrayList<Integer>();
	ArrayList<Integer> visitedWards = new ArrayList<Integer>(); // TODO record visited LSOAs
	ArrayList<Polygon> polys = new ArrayList<Polygon>();
	ArrayList<String> csvData = new ArrayList<String>();

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
			ShapeFileImporter.read(roadsFile.toURI().toURL(), roadLayer);
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

			//createNetwork();

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
			roadLayer.setMBR(MBR);
			flood3.setMBR(MBR);
			flood2.setMBR(MBR);
			
			baseLayer.setMBR(MBR);

			// clean up the road network
			
			System.out.print("Cleaning the road network...");
			
			roadNetwork = NetworkUtilities.multipartNetworkCleanup(roadLayer, roadNodes, resolution, fa, random, 0);
			roadNodes = roadNetwork.getAllNodes();
			//testNetworkForIssues(roadNetwork);
			
			// set up roads as being "open" and assemble the list of potential terminii
			roadLayer = new GeomVectorField(grid_width, grid_height);
			for(Object o: roadNodes){
				GeoNode n = (GeoNode) o;
				networkLayer.addGeometry(n);
				
				boolean potential_terminus = false;
				
				// check all roads out of the nodes
				for(Object ed: roadNetwork.getEdgesOut(n)){
					
					// set it as being (initially, at least) "open"
					ListEdge edge = (ListEdge) ed;
					((MasonGeometry)edge.info).addStringAttribute("open", "OPEN");
					networkEdgeLayer.addGeometry( (MasonGeometry) edge.info);
					roadLayer.addGeometry((MasonGeometry) edge.info);
					((MasonGeometry)edge.info).addAttribute("ListEdge", edge);
					
					String type = ((MasonGeometry)edge.info).getStringAttribute("TYPE");
					if(type.equals("motorway") || type.equals("primary") || type.equals("trunk"))
						potential_terminus = true;
				}
				
				// check to see if it's a terminus
				if(potential_terminus && !MBR.contains(n.geometry.getCoordinate()) && roadNetwork.getEdges(n, null).size() == 1){
					terminus_points.add(n);
				}

			}
			
			// reset MBRS in case it got messed up during all the manipulation
			baseLayer.setMBR(MBR);
			roadLayer.setMBR(MBR);
			flood3.setMBR(MBR);
			flood2.setMBR(MBR);
			
			System.out.println("Done cleaning the road network!");
			
			/////////////////////
			///////// Clean up roads for Agents to use ///////////
			/////////////////////
						
			Network majorRoads = extractMajorRoads();
			testNetworkForIssues(majorRoads);

			// assemble list of secondary versus local roads
			ArrayList <Edge> myEdges = new ArrayList <Edge> ();
			GeomVectorField secondaryRoadsLayer = new GeomVectorField(grid_width, grid_height);
			GeomVectorField localRoadsLayer = new GeomVectorField(grid_width, grid_height);
			for(Object o: majorRoads.allNodes){
				
				majorRoadNodesLayer.addGeometry((GeoNode)o);
				
				for(Object e: roadNetwork.getEdges(o, null)){
					Edge ed = (Edge) e;
					
					myEdges.add(ed);
										
					String type = ((MasonGeometry)ed.getInfo()).getStringAttribute("class");
					if(type.equals("secondary"))
							secondaryRoadsLayer.addGeometry((MasonGeometry) ed.getInfo());
					else if(type.equals("local"))
							localRoadsLayer.addGeometry((MasonGeometry) ed.getInfo());					
				}
			}

			System.gc();
			
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
	 * Connect the GeoNode to the given subnetwork using the complete road network
	 * 
	 * @param n - the target node
	 * @param subNetwork - the existing subnetwork
	 */
	void connectToMajorNetwork(GeoNode n, Network subNetwork) {

		try {
			Bag subNetNodes;			
			subNetNodes = (Bag) subNetwork.allNodes.clone();
			
			// find a path using the whole set of roads in the environment 
			AStar pathfinder = new AStar();
			ArrayList <Edge> edges = pathfinder.astarPath(n, new ArrayList <GeoNode> (subNetNodes), roadNetwork);
			
			if(edges == null) return; // maybe no such path exists!

			//  otherwise, add the edges into the subnetwork
			for(Edge e: edges){
				GeoNode a = (GeoNode) e.getFrom(), b = (GeoNode) e.getTo();
				if(!subNetwork.nodeExists(a) || !subNetwork.nodeExists(b))
					subNetwork.addEdge(a, b, e.info);
			}

		} catch (CloneNotSupportedException e1) {
			e1.printStackTrace();
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
	 * Make sure the network doesn't have any problems
	 * 
	 * @param n - the network to be tested
	 */
	static void testNetworkForIssues(Network n){
		System.out.println("testing");
		for(Object o: n.allNodes){
			GeoNode node = (GeoNode) o;
			for(Object p: n.getEdgesOut(node)){
				sim.field.network.Edge e = (sim.field.network.Edge) p;
				LineString ls = (LineString)((MasonGeometry)e.info).geometry;
				Coordinate c1 = ls.getCoordinateN(0);
				Coordinate c2 = ls.getCoordinateN(ls.getNumPoints()-1);
				GeoNode g1 = (GeoNode) e.getFrom();
				GeoNode g2 = (GeoNode) e.getTo();
				if(c1.distance(g1.geometry.getCoordinate()) > 1)
					System.out.println("found you");
				if(c2.distance(g2.geometry.getCoordinate()) > 1)
					System.out.println("found you");
			}
		}
	}
	
	/**
	 * //////////////////////// Create Road Network //////////////////////////////
	 * Create the road network the agents will traverse
	 */
	/*
	private void createNetwork() {
		System.out.println("Creating road network..." + roadLayer);
		System.out.println();
		network.createFromGeomField(roadLayer);

		for (Object o : network.getEdges()) {
			GeomPlanarGraphEdge e = (GeomPlanarGraphEdge) o;

			idsToEdges.put(e.getIntegerAttribute("ROAD_ID_1").intValue(), e);

			e.setData(new ArrayList<agents.Agent>());
		}

		addIntersectionNodes(network.nodeIterator(), junctions);
	}
	*/

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

			BufferedReader agentData = new BufferedReader(new InputStreamReader(fstream));
			String s;

			// get rid of the header
			agentData.readLine();
			// read in all data
			while ((s = agentData.readLine()) != null) {
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

					MasonGeometry newGeometry = new MasonGeometry(a.getGeometry());
					//MasonGeometry newGeometry = a.getGeometry();
					newGeometry.isMovable = true;
					agentsLayer.addGeometry(newGeometry);
					this.agentList.add(a);
					schedule.scheduleOnce(a);
				}
			}

			agentData.close();
			System.out.println();
			System.out.println("All Agents added successfully!");
		} catch (Exception e) {
			System.out.println("ERROR: issue with population file: ");
			e.printStackTrace();
		}
	}
	
	/**
	 * Coordinate reader helper function
	 * @param s
	 * @return
	 */
	Coordinate readCoordinateFromFile(String s){
		if(s.equals("")) 
			return null;
		
		String [] bits = s.split(" ");
		Double x = Double.parseDouble( bits[1].substring(1) );
		Double y = Double.parseDouble(bits[2].substring(0, bits[2].length() - 2));
		return new Coordinate(x,y);
	}
	
	/**
	 * Extract the major roads from the road network
	 * @return a connected network of major roads
	 */
	public Network extractMajorRoads(){
		Network majorRoads = new Network();
		
		// go through all nodes
		for(Object o: roadNetwork.getAllNodes()){
		
			GeoNode n = (GeoNode) o;
			
			// go through all edges
			for(Object p: roadNetwork.getEdgesOut(n)){
				
				sim.field.network.Edge e = (sim.field.network.Edge) p;
				String type = ((MasonGeometry)e.info).getStringAttribute("class");
				
				// save major roads
				if(type.equals("major"))
						majorRoads.addEdge(e.from(), e.to(), e.info);
			}
		}
		
		// merge the major roads into a connected component
		NetworkUtilities.attachUnconnectedComponents(majorRoads, roadNetwork);
		
		return majorRoads;
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
		Bag objects = roadLayer.getObjectsWithinDistance(fa.createPoint(c), resolution);
		if (objects == null || roadLayer.getGeometries().size() <= 0)
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
		Bag objects = roadLayer.getObjectsWithinDistance(fa.createPoint(c), resolution);
		if (objects == null || roadLayer.getGeometries().size() <= 0)
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
		Bag objects = roadLayer.getObjectsWithinDistance(fa.createPoint(c), resolution);
		if (objects == null || roadLayer.getGeometries().size() <= 0)
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
	
	// reset the agent layer's MBR
		public void resetLayers(){
			MBR = baseLayer.getMBR();
			MBR.init(501370, 521370, 4292000, 4312000);
			this.agentsLayer.setMBR(MBR);
			this.roadLayer.setMBR(MBR);
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
