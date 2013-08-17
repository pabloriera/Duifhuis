#pragma once

#include "ofMain.h"
#include <iostream>
#include <vector>
#include "../src/cochlea.h"


using namespace std;
//using namespace arma;

class testApp : public ofBaseApp{
	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y);
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);

		double *Y_t, *V_t;
        int n_t,n;
        int n_osc;
        float fs;
        Cochlea_t cochlea;

        double* stimulus;

        float A;
        float dL;
        float L;
        float zscale;

        //Stimulus
        float f0, sig, onset_dur;
        float t;


};
