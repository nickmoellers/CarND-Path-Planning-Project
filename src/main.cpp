#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"

#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}
/*
vector<double> make_ptss( vector<double> ptsx, vector<double> ptsy, double theta,
  const vector<double> &maps_x, const vector<double> &maps_y) {
  
  double x,y;
  double prev_x, prev_y;
  vector<double> sd;
  vector<double> ptss;

  for( int i = 0 ; i < ptsx.size() ; i ++ ) {
    if( i > 0 ) theta = atan2(y-prev_y, x-prev_x); 

    x = ptsx[i];
    y = ptsy[i];

    sd = getFrenet( x,  y,  theta, maps_x, maps_y );
    double s = sd[0];
    ptss.push_back(s);

    prev_x = x;
    prev_y = y;
  }
  return ptss;
}*/

int dToLane( double d ) {
  int lane = floor(d)/4;
  return lane;
}

double laneToD( int lane ) {
  double d = 2+4*lane;
  return d;
}

vector<double> check_lane( int car_lane, int check_direction, vector<double> my_car, vector<vector<double>> sensor_fusion ) {

  //cout << "car_lane: " << car_lane << endl;
  
  /*double min_time_to_collision = 100000.0; // large number
  double min_distance_to_collision = 100000.0; //large number
  double collision_id = -1;*/
  double next_car_id = -1;
  double next_car_delta_s = check_direction * 10000000; //large numbers
  double next_car_delta_vs = check_direction * 1000000;
  double next_car_s = -1;
  double next_car_vs = -1;

  //double car_id = my_car[0]; //-1;
  double car_x =  my_car[1];
  double car_y =  my_car[2];
  double car_vx = my_car[3];
  double car_vy = my_car[4];
  double car_s =  my_car[5];
  double car_d =  my_car[6];
  //assume vd is 0 so vs is the hypotenuse of vx and vy and no lane change.
  double car_vs_mph = distance(0,0,car_vx,car_vy); // mph
  double car_vs = car_vs_mph / 2.24; //m/s
  //double car_lane = dToLane(d);

  cout << "my_car:" << endl;
  // cout << "\tx:\t" << car_x << endl;
  // cout << "\ty:\t" << car_y << endl;    
  // cout << "\tvx:\t" << car_vx << endl;
  // cout << "\tvy:\t" << car_vy << endl;
  cout << "\ts:\t" << car_s << endl;
  // cout << "\td:\t" << car_d << endl;
  cout << "\tvs:\t" << car_vs << endl;
  // cout << "\tlane:\t" << car_lane << endl << endl;
  // cout << "\tcalculated_lane:\t" << dToLane(car_d) << endl << endl;

  double id, x, y, vx, vy, s, d, vs, lane;
  double delta_s, delta_vs;
  for( vector<double> sensor_fusion_i : sensor_fusion ) {
    id = sensor_fusion_i[0];
    x = sensor_fusion_i[1];
    y = sensor_fusion_i[2];
    vx = sensor_fusion_i[3];
    vy = sensor_fusion_i[4];
    s = sensor_fusion_i[5];
    d = sensor_fusion_i[6];
    //assume vd is 0 so vs is the hypotenuse of vx and vy and no lane change.
    vs = distance(0,0,vx,vy);
    lane = dToLane(d);

    delta_s = 0;
    delta_vs = 0;

    if( car_lane == lane ) {

      

      delta_s = s-car_s; //positive for car in front
      delta_vs = car_vs-vs; //postiive if I'm going faster
      //multiple by check direction (which is positive if checking for car in front)
      double target_direction = check_direction * delta_s;

      if( target_direction > 0 ) {
        if( abs(delta_s) < abs(next_car_delta_s) ) {
          next_car_id = id;
          next_car_delta_s = delta_s;
          next_car_delta_vs = delta_vs;
          next_car_s = car_s;
          next_car_vs = car_vs;
          cout << "*";
        }
        cout << "car_" << id << ":" << endl;
        // cout << "\tx:\t" << x << endl;
        // cout << "\ty:\t" << y << endl;
        // cout << "\tvx:\t" << vx << endl;
        // cout << "\tvy:\t" << vy << endl;
        cout << "\ts:\t" << s << endl;
        // cout << "\td:\t" << d << endl;
        cout << "\tvs:\t" << vs << endl;
        //cout << "\tlane:\t" << lane << endl << endl;
        //cout << "\tcar_lane:\t" << car_lane << endl;


        // cout << "check_direction: " << check_direction << endl;
        // cout << "delta_s: " << delta_s << endl;
        // cout << "delta_vs: " << delta_vs << endl;
        // cout << "target_direction: " << target_direction << endl;

        
      }
    }
  }

  vector<double> next_car; 
  next_car.push_back(next_car_id); //0
  //collision.push_back(min_time_to_collision); //1
  next_car.push_back(next_car_delta_s);
  next_car.push_back(next_car_delta_vs);
  next_car.push_back(next_car_s);
  next_car.push_back(next_car_vs);
  //collision.push_back(car_s);
  //collision.push_back(car_vs);
  //collision.push_back(s);
  //collision.push_back(vs);
  
  return next_car;
}

int change_lanes( int this_lane , double next_car_delta_s, vector<double> car, vector<vector<double>> sensor_fusion ) {
  cout << "Should I change into lane " << this_lane << "?" << endl;
  if( this_lane < 0 || this_lane > 2 ) {
    cout << "no such lane!" << endl;
    return -1;
  }
  // int change_lane = lane-1;
  // if( change_lane >= 0 && change_lane <= 2 ) {
  //   cout << "Checking lane: " << change_lane << endl;
  vector<double> change_lane_car_fwd = check_lane( this_lane, 1, car, sensor_fusion );
  double change_lane_car_fwd_delta_s = change_lane_car_fwd[1];

  vector<double> change_lane_car_back = check_lane( this_lane, -1, car, sensor_fusion );
  double change_lane_car_back_delta_s = change_lane_car_back[1];
  
  cout << "There's " << change_lane_car_fwd_delta_s << " meters free in front of me and" << endl;
  cout << change_lane_car_back_delta_s << " meters free behind me." << endl;
  if( change_lane_car_fwd_delta_s > next_car_delta_s + 10 && 
      change_lane_car_fwd_delta_s > 50 &&
      change_lane_car_back_delta_s < 10 ) { //don't cut someone off!
    cout << "Yes, change lanes!" << endl;
    return this_lane;
  } else {
    cout<< "Nope, lane not clear." << endl;
    return -1;
  }
  //}

  
}

double speed_limit_mph = 50; // mph
  double speed_limit = speed_limit_mph/2.24; //m/s
  double target_velocity = speed_limit;
  double dt = 0.02; // seconds
  int lane = 1;

int main() {
  uWS::Hub h;

  //int lane = 1;

  


  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          ////////////////////////////////////////////////////////////
         
          // cout << "previous_path(x, y): " << previous_path_x.size() << endl;
          // for( int i = 0 ; i < previous_path_x.size() ; i ++ ) {
          //   cout << i << ":\t(" << previous_path_x[i] << ",\t" << previous_path_y[i] << ")" << endl;
          // }

          //int lane = 1;
          //double speed_limit = 50.0;

          double ref_x, ref_y;
          double ref_yaw = deg2rad(car_yaw);
          double prev_x, prev_y;
          double reference_velocity;

          //
          if( previous_path_x.size()<=1 ) {
            ref_x = car_x;
            ref_y = car_y;

            prev_x = car_x - cos(car_yaw);
            prev_y = car_y - sin(car_yaw);

            reference_velocity = car_speed;
          } else {
            ref_x = previous_path_x[previous_path_x.size()-1];
            ref_y = previous_path_y[previous_path_x.size()-1];
            
            prev_x = previous_path_x[previous_path_x.size()-2];
            prev_y = previous_path_y[previous_path_y.size()-2];

            ref_yaw = atan2(ref_y-prev_y, ref_x-prev_x);

            reference_velocity = distance(prev_x, prev_y, ref_x, ref_y)/dt;
          

            ///CHECK LANE HERE
            /*vector <double> ptss;
            ptss = make_ptss( previous_path_x, previous_path_y, deg2rad(car_yaw), map_waypoints_x, map_waypoints_y );
            cout << "ptss: " << ptss.size() << endl;
            for( int i = 0 ; i < ptss.size() ; i ++ ) {
              cout << i << ":\t(" << ptss[i] << " )" << endl;
            }*/

            vector<double> car;
            car.push_back(-1.0); //id;
            car.push_back(car_x); //x;
            car.push_back(car_y); //y
            double vx = car_speed * cos(car_yaw);
            double vy = car_speed * sin(car_yaw);
            //yeah, So I could just lie and pass in car_speed and 0, but why not?
            car.push_back(vx); //vx
            car.push_back(vy); //vy
            car.push_back(car_s); //s
            car.push_back(car_d); //d

            vector<double> next_car = check_lane( lane, 1, car, sensor_fusion );
            int next_car_id = next_car[0];
            double next_car_delta_s = next_car[1];
            double next_car_delta_vs = next_car[2];
            double next_car_s = next_car[3];
            double next_car_vs = next_car[4];

            if( next_car_id > 0) {
              cout << "the next vehicle in my lane, car # " << next_car_id;
              cout << " is currently " << next_car_delta_s << " meters away." << endl;
              cout << "I am approaching it at a relative speed of " << next_car_delta_vs << endl;
              double next_car_collide_time = next_car_delta_s / next_car_delta_vs;
              cout << "I will hit it in " << next_car_collide_time << " seconds. " << endl;

              if( next_car_delta_s < 30 ) { //If I am close to a car
                cout << "CLOSE TO CAR" << endl;
                if( next_car_delta_vs > -1 ) { //If I am going faster than that car 
                  cout << "APPROACHING CLOSE CAR" << endl;
                  //target_velocity-=next_car_delta_vs; //next_car_delta_vs;
                  target_velocity = next_car_vs-5;//max( next_car_vs - 5, target_velocity );
                  
                } //else { // If I am going slower than that car
              } else if ( next_car_delta_s > 90 ) {
                cout << "I am greater than 90 meters behind the car " << endl;
                target_velocity = speed_limit-1;
              } else  {
                // maybe change lanes?
                // check left and right
                int new_lane = change_lanes( lane-1 , next_car_delta_s, car, sensor_fusion ) ;
                if( new_lane < 0  ) {
                  new_lane = change_lanes( lane+1 , next_car_delta_s, car, sensor_fusion ) ;
                }


                if( new_lane < 0  ) {
                  //can't change lanes!
                  /*// follow vehicle at 60 meters back
                target_velocity+= 0.5;
                target_velocity = max( next_car_vs - 2, target_velocity );*/
                } else {
                  lane = new_lane;
                }
                

                
             
              }
            } else {
              cout << "lane comletely empty in sensor horizon" << endl;
              target_velocity = speed_limit-1;
            }
            cout << "target_velocity: " << target_velocity << endl;
/*
              cout << "collision with " << collision_id << " in t = " << time_to_collision << "seconds" << endl;
              if( abs( time_to_collision ) < 2 ) { //less than 2 seconds is tailgating
                cout << "collision soon!" << endl;
                double collision_car_vx = sensor_fusion[collision_id][3];
                double collision_car_vy = sensor_fusion[collision_id][4];
                double vs = distance(0,0,vx,vy);
                cout << "reducing speed to: " << vs << endl;
                double target_velocity_mph = vs;
                target_velocity = target_velocity_mph/2.24-2;
              } else if ( abs( time_to_collision ) < 10 )  {
                cout << "pacing car..." << endl;
                cout << "increasing speed slightly..." << endl;
                target_velocity+=0.5;
              } else {
                cout << "collision vehicle very far!" << endl;
                target_velocity = speed_limit -1;
              }
              //if( abs( time_to_collision ) )
            } else {
              cout << "no collsion vehicles!" << endl;
              target_velocity = speed_limit -1;
            }
           */
            //no harm in checking...
            if( target_velocity > speed_limit -1 ) {
              target_velocity = speed_limit - 1;
            }

          

          }

          vector<double> ptsx;
          vector<double> ptsy;

          ptsx.push_back(prev_x);
          ptsy.push_back(prev_y);

          ptsx.push_back(ref_x);
          ptsy.push_back(ref_y);

          // cout << "car:\t(" << car_x << ",\t" << car_y << "), yaw =\t" << car_yaw << endl;
          // cout << "ref:\t(" << ref_x << ",\t" << ref_y << "), yaw =\t" << ref_yaw << endl;
          // cout << "prev:\t(" << prev_x << ",\t" << prev_y << ")" << endl;

          // cout << "pts(x, y): " << ptsx.size() << endl;
          // for( int i = 0 ; i < ptsx.size() ; i ++ ) {
          //   cout << i << ":\t(" << ptsx[i] << ",\t" << ptsy[i] << ")" << endl;
          // }
            
          //add points to spline at 30, 60, 90 meters out
          for( int i = 50; i<=100 ; i+=25 ) {
            vector<double> next_xy = getXY(car_s + i, laneToD(lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            ptsx.push_back(next_xy[0]);
            ptsy.push_back(next_xy[1]);
          }

          // cout << "pts(x, y): " << ptsx.size() << endl;
          // for( int i = 0 ; i < ptsx.size() ; i ++ ) {
          //   cout << i << ":\t(" << ptsx[i] << ",\t" << ptsy[i] << ")" << endl;
          // }

          //convert to car's reference frame
          for( int i = 0; i < ptsx.size(); i++ ) {
            double shift_x = ptsx[i]-ref_x;
            double shift_y = ptsy[i]-ref_y;
            ptsx[i] = shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw);
            ptsy[i] = shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw);
          }

          // cout << "pts(x, y): " << ptsx.size() << endl;
          // for( int i = 0 ; i < ptsx.size() ; i ++ ) {
          //   cout << i << ":\t(" << ptsx[i] << ",\t" << ptsy[i] << ")" << endl;
          // }

          tk::spline s;
          s.set_points(ptsx,ptsy);

          for( int i = 0 ; i < previous_path_x.size() ; i++ ) {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          // cout << "next_vals_from_prev_path(x, y): " << next_x_vals.size() << endl;
          // for( int i = 0 ; i < next_x_vals.size() ; i ++ ) {
          //   cout << i << ":\t(" << next_x_vals[i] << ",\t" << next_y_vals[i] << ")" << endl;
          // }

        

          

          vector<double> x_spline_vals;
          vector<double> y_spline_vals;

          double x_spline;
          double y_spline;
          double prev_x_spline = 0;

          //cout << "speed_limit: " << speed_limit << endl;
          //cout << "target_velocity: " << target_velocity << endl;
          
          while( next_x_vals.size() < 50 ) {


            //cout << "\treference_velocity: " << reference_velocity << endl;

            if( reference_velocity < target_velocity -1 &&
                reference_velocity < speed_limit -1 )
                reference_velocity+=0.19;
            if( reference_velocity > target_velocity + 1 )
                reference_velocity-=0.19;

            //double reference_velocity_m_s = reference_velocity/2.24; //m/s
            //double dt = 0.02; // seconds
            double target_x = 30; // m
            double target_y = s(target_x); //m
            double target_dist = distance(0, 0, target_x, target_y); //m
            double N = target_dist / (0.02*reference_velocity); //unitless
            double target_distance_per_timestep = target_x/N; // m
            //cout << "deltaX/t = " << target_distance_per_timestep << endl;

            x_spline = prev_x_spline + target_distance_per_timestep; 
            y_spline = s(x_spline);

            prev_x_spline = x_spline;
            
            double x_calc = ( x_spline*cos(ref_yaw) - y_spline*sin(ref_yaw) );
            double y_calc = ( x_spline*sin(ref_yaw) + y_spline*cos(ref_yaw) );
            
            double x_point = x_calc + ref_x;
            double y_point = y_calc + ref_y;

            x_spline_vals.push_back(x_point);
            y_spline_vals.push_back(y_point);

            next_x_vals.push_back(x_point);
            next_y_vals.push_back(y_point);

          }

          // cout << "next_vals(x, y): " << next_x_vals.size() << endl;
          // for( int i = 0 ; i < next_x_vals.size() ; i ++ ) {
          //   cout << i << ":\t(" << next_x_vals[i] << ",\t" << next_y_vals[i] << ")" << endl;
          // }

          ////////////////////////////////////////////////////////////

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          //cout << "MESSAGE SENDING" << endl;
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          //cout << "MESSAGE SENT" << endl;
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen("0.0.0.0",port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
