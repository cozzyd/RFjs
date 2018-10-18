"use strict"; 

/** 
 *  Various things for working with RF TGraphs
 *   Sorry, I dont' really know javascript, so this is probably all terrible. 
 *
 */ 

var RF = {}; 
RF.ffts = {}; 

/** Obtain an FFT of the right size 
 *
 * Then you an call forward or inverse on the returned object. 
 *
 **/ 
RF.getFFT =function (size) 
{
  if (this.ffts[size] === undefined)
  {
    this.ffts[size] = new FFTR(size); 
  }

  return this.ffts[size]; 
};


/** performs FFT */ 

RF.doFFT = function (y) 
{
    var fft = RF.getFFT(y.length); 
    return fft.forward(y).slice(0); 
};


/** Performs inverse FFT */ 
RF.doInvFFT = function(Y, Nt = 0) 
{
  var fft = RF.getFFT(Nt == 0 ? Y.length : Nt); 
  return fft.inverse(Y).slice(0); 
};


/** Returns the maximum time and value */ 
RF.getMaximumTimeAndValue = function(g, unsigned=true) 
{
  var max = 0; 
  var max_t = -1;

  for (var i = 0; i < g.fNpoints; i++) 
  {
    var val = g.fY[i]; 
    if (unsigned && val < 0) val = -val; 
    
    if (val > Math.abs(max)) 
    {
      max = g.fY[i]; 
      max_t = g.fX[i]; 
    }
  }

  return [max_t, max]; 
}

/** Generates a power spectrum of the given TGraph, returning another TGraph.
 * If you already have the FFT available, you can pass it as Y
 * */ 

RF.makePowerSpectrum = function(g, Y = null) 
{


  var N = g.fX.length;
  if (Y == null) 
  {
    var fft = RF.getFFT(N); 
    Y = fft.forward(g.fY).slice(0); 
  }
  var dt = g.fX[1] - g.fX[0]; 
  var df = 1./(N * dt); 
  var f = []; 
  var P = []; 

  for (var i = 0; i < N/2+1; i++)
  {
    f.push(i*df); 
    var p = (Y[2*i]*Y[2*i] + Y[2*i+1]*Y[2*i+1]) / N; 
    if (i > 0 || i < N/2) p*=2; 
    P.push(10*Math.log10(p)); 
  }

  var G = JSROOT.CreateTGraph(N/2+1, f,P); 
  G.fName = g.fName + "_power"; 
  G.fTitle = g.fTitle.substring(g.fTitle.indexOf(',')+1); 
 
  G.InvertBit(JSROOT.BIT(18)); //make it uneditable

  return G;
}




/** upsamples (in place) a TGraph in the fourier domain 
 *
 *  g: The TGraph
 *  factor: the upsample factor (>1) 
 *  Y: If you already have the fourier trans
 *
 *  returns the fourier transform of the new graph (in case you want to do more with it) 
 *
 * */ 
RF.upsample = function (g, factor=2, Y = null) 
{
  if (factor <=1) return RF.doFFT(g.fY); 

  var N = g.fNpoints; 
  var t0 = g.fX[0]; 
  var dt = g.fX[1] - g.fX[0]; 
  if (Y == null) 
  {
    Y = RF.doFFT(g.fY); 
  }

  var newY = new Float32Array( 2*(g.fNpoints * factor / 2 + 1)); 

  for (var i = 0; i < Y.length; i++) 
  {
    newY[i] = Y[i]/N; 
  }

  g.fNpoints*= factor; 
  g.fY =RF.doInvFFT(newY); 

  for (var i = 0; i < g.fNpoints; i++ )
  {
    g.fX[i] = dt / factor * i + t0 ;
  }
  
  return newY;
}


/** Returns the hilbert transform of g */ 
RF.hilbertTransform = function(g, Y = null) 
{

  if (Y == null) Y = RF.doFFT(g.fY); 
  var Yp = new Float32Array(Y.length); 

  var N = g.fNpoints; 
  for (var i = 1; i < g.fNpoints/2; i++) 
  {

    Yp[2*i] = -Y[2*i+1] / N; 
    Yp[2*i+1] = Y[2*i] / N; 
  }

  var yp = RF.doInvFFT(Yp); 
  var gh= JSROOT.CreateTGraph(g.fNpoints, g.fX, yp); 
  gh.fName = g.fName + "_hilbert"; 
  gh.fTitle = "Hilbert Transform of " + g.fTitle; 
  return gh; 
}


/** Creates the hilbert envelope of a graph */ 

RF.hilbertEnvelope = function(g, Y = null, gh = null, H = null) 
{
  if (Y == null) Y = RF.doFFT(g.fY); 

  if (gh == null) 
  {
    gh = RF.hilbertTransform(g,Y); 
  }

  if (H == null) 
  {
    H = JSROOT.CreateTGraph(g.fNpoints, g.fX, g.fY); 
  }
  else
  {
    H.fNpoints = g.fNpoints; 
    H.fX = g.fX; 
    H.fY = []; 
  }

  H.fName = g.fName+"_envelope"; 
  H.fTitle = "Hilbert envelope of " + g.fTitle;
  for (var i = 0; i < g.fNpoints; i++)
  {
    H.fY[i] = Math.sqrt( g.fY[i] * g.fY[i] + gh.fY[i] * gh.fY[i]); 
  }

  H.InvertBit(JSROOT.BIT(18)); //make it uneditable
  return H ; 
}



RF.getGraphMean = function(g) 
{
  var sum = 0; 
  for (var i = 0; i < g.fNpoints; i++) sum += g.fY[i]; 
  return sum / g.fNpoints; 
}

RF.getMean = function(y) 
{
  var sum = 0; 
  for (var i = 0; i < y.length; i++) sum += y[i]; 
  return sum / y.length; 
}

RF.getRMS = function(g) 
{
  var sum2 = 0; 
  var sum = 0; 
  var N = g.fNpoints;

  for (var i = 0; i < g.fNpoints; i++) 
  {
    sum += g.fY[i];
    sum2 += g.fY[i] * g.fY[i];
  }
  var mean = sum/N;

  return Math.sqrt(sum2/N-mean*mean);
}




RF.rectifyGraph = function(g) 
{
  var mean = RF.getGraphMean(g); 
  for (var i = 0; i < g.fNpoints; i++) g.fY[i] -= mean; 

}

RF.rectify= function(y) 
{
  var mean = RF.getMean(y); 
  for (var i = 0; i < y.length; i++) y[i]-=mean; 

}





RF.range = function(start, N, step = 1) 
{
  var ans = new Float32Array(N); 
  for (var i = 0; i < N; i++) 
  {
    ans[i] = start+step*i; 
  }

  return ans; 
}


/* Calculates the cross-correlation of two graphs, returning the cross-correlation the time domain 
 * Asumes they have the same sampling and start position 
 * */ 

RF.crossCorrelation = function ( g1, g2, pad=true, upsample = 4, scale = null, cutoff = 0) 
{
  
  if (scale == null) scale = RF.getRMS(g1) * RF.getRMS(g2) ;


  /** most likely we need to pad by a factor of 2 */ 

  var u1 = RF.getGraphMean(g1);
  var u2 = RF.getGraphMean(g2);
  var N = Math.max(g1.fNpoints, g2.fNpoints); 
  var dt = g1.fX[1] - g1.fX[0]; 
  if (pad) N*=2; 
  
  scale *= N/2 * N; 
  var df = 1./(N*dt); 

  var y1 = new Float32Array(N); 
  var y2 = new Float32Array(N); 

  for (var i = 0; i < g1.fNpoints;i++) y1[i] = g1.fY[i]-u1; 
  for (var i = 0; i < g2.fNpoints;i++) y2[i] = g2.fY[i]-u2; 


  var Y1 = RF.doFFT(y1); 
  var Y2 = RF.doFFT(y2); 

  var Y = new Float32Array(upsample * Y1.length); 

  for (var i = 0; i < N/2+1;i++) 
  {
    if (cutoff>0 && i*df > cutoff) break; 

    var re1 = Y1[2*i];
    var re2 = Y2[2*i];
    var im1 = Y1[2*i+1];
    var im2 = Y2[2*i+1]; 

    Y[2*i] = (re1*re2 + im1*im2)/scale;
    Y[2*i+1] = (im1*re2-re1*im2)/scale;
  }

  var y = RF.doInvFFT(Y); 
  var yrotated = new Float32Array (y.length); 
  for (var i = 0; i < y.length; i++) 
  {
    yrotated[i] = y[(i+y.length/2) % y.length];
  }

  N = y.length; 
  dt = dt/upsample; 
  var x = RF.range(-N/2*dt, N, dt); 
  var g = JSROOT.CreateTGraph(y.length, x, yrotated); 

  g.fYitle = "Correlation of" + g1.fTitle + " with " + g2.fTitle; 

  return g; 
}


/** Interpolates a graph at a point, assuming it's evenly spaced */ 
RF.evalEven = function (g, t) 
{

  var t0 = g.fX[0]; 
  var dt = g.fX[1] - g.fX[0]; 
  var l = Math.floor((t-t0)/dt);
  var u = l+1; 
  if (l < 0) return 0; 
  if (l >= g.fNpoints) return 0; 

  var f = (t - (t0+dt*l))/dt;  
  return f * g.fY[u] + (1-f) * g.fY[l]; 
}



RF.Mapper = function (nants, computeDeltaTs, usePair = function(i,j) { return true; }, canUse = function(i, x, y) { return true; } ) 
{
  var m = {}; 
  m.nant = nants; 
  m.deltaTs = computeDeltaTs; 
  m.usePair = usePair; 
  m.canUse = canUse; 
  return m; 
}




RF.Antenna = function(x,y,z, dx,dy,dz, max_phi,max_theta) 
{
  var ant = {};
  ant.pos = [x,y,z];
  ant.bore = [dx,dy,dz]; 
  ant.phi_width = max_phi;
  ant.theta_width = max_theta; 
  return ant; 
}


RF.dotProduct = function(x,y, Nmax = 0) 
{
  if (Nmax == 0) Nmax = x.length; 
  var sum = 0; 
  for (var i = 0; i < Nmax; i++) sum+= x[i] * y[i]; 
  return sum; 
}

RF.crossProduct = function(u,v) 
{

  return [ u[1]*v[2] - u[2]*v[1] ,
           u[2]*v[0] - u[0]*v[2] , 
           u[0] *v[1] - u[1]*v[0] ]; 
}


RF.magnitude = function(x, Nmax) 
{
  return Math.sqrt(RF.dotProduct(x,x,Nmax)); 
}

/* Returns the angle between two vectors, in radians */ 
RF.angleBetween = function(x,y) 
{
  return Math.acos( RF.dotProduct(x,y) / (RF.magnitude(x) * RF.magnitude(y))); 
}


/* From  https://stackoverflow.com/questions/966225/how-can-i-create-a-two-dimensional-array-in-javascript */ 
RF.createArray = function(length) {
    var arr = new Array(length || 0),
        i = length;
    if (arguments.length > 1) {
        var args = Array.prototype.slice.call(arguments, 1);
        while(i--) arr[length-1 - i] = RF.createArray.apply(this, args);
    }
    return arr;
}



RF.wrap = function(x, period = 360, center = 0)
{
  return x - period * Math.floor((x-center+period/2)/period); 
}

RF.AngleMapper = function ( ants,  c = 0.3, phi_0 =0, theta_0 = 0 )
{
  var n = ants.length; 
  var can_use = RF.createArray(n,n); 
  var min_phi =  RF.createArray(n); 
  var max_phi = RF.createArray(n); 
  var min_theta = RF.createArray(n); 
  var max_theta = RF.createArray(n);

  for (var i= 0; i< n; i++)
  {
    var bore_phi = 180./Math.PI* Math.atan2(ants[i].bore[1], ants[i].bore[0]) + phi_0; //check the sign
    var bore_theta = 180./Math.PI* Math.asin(ants[i].bore[2]) + theta_0; //check the sign

    if (ants[i].phi_width >=180) 
    {
      min_phi[i] =-180; 
      max_phi[i] = 180 
    }
    else
    {
      min_phi[i] = RF.wrap(bore_phi - ants[i].phi_width,360,0); 
      max_phi[i] = RF.wrap(bore_phi + ants[i].phi_width,360,0); 
    }

    if (ants[i].theta_width >= 90)
    {

      min_theta[i] =-90; 
      max_theta[i] = 90;
 
    }
    else
    {
      min_theta[i] = RF.wrap(bore_theta - ants[i].theta_width,180,0); 
      max_theta[i] = RF.wrap(bore_theta + ants[i].theta_width,180,0); 
    }

    can_use[i][i] = false; 

    for (var j = i+1; j < n; j++) 
    {
      var bore_xy_i = [ants[i].bore[0], ants[i].bore[1], 0] ; 
      var bore_xy_j = [ants[j].bore[0], ants[j].bore[1], 0] ; 
      var bore_rz_i = [RF.magnitude(ants[i].bore,2),0, ants[i].bore[2]] ; 
      var bore_rz_j = [RF.magnitude(ants[j].bore,2),0, ants[j].bore[2]] ; 
      var boresight_phi_angle =  RF.angleBetween(bore_xy_i, bore_xy_j) * 180/Math.PI; 
      var boresight_theta_angle =  RF.angleBetween(bore_rz_i, bore_rz_j) * 180/Math.PI; 

      can_use[i][j] = boresight_phi_angle < ants[i].phi_width && boresight_phi_angle < ants[j].phi_width && boresight_theta_angle < ants[i].theta_width && boresight_theta_angle < ants[j].theta_width; 
      can_use[j][i] = can_use[i][j]; 
//      console.log(i,j,can_use[i][j]); 
    }
  }
   

  return RF.Mapper( ants.length , 

    function (i,j, phi_deg, theta_deg)  //deltaTs
    {
      var phi = Math.PI / 180. * (phi_deg - phi_0); 
      var theta = Math.PI / 180. * (theta_deg - theta_0); 
      var dir = [ Math.cos(phi) * Math.cos(theta), Math.sin(phi) * Math.cos(theta), Math.sin(theta) ]; 
      var diff = [ ants[j].pos[0] - ants[i].pos[0], ants[j].pos[1] - ants[i].pos[1], ants[j].pos[2] - ants[i].pos[2] ];
      return RF.dotProduct(diff,dir) / c; 
    }, 

    function( i, j)  //usePair
    {
      return can_use[i][j]; 

    }, 

    function(i,phi_deg,theta_deg)  //canUse
    {

       if (theta_deg > max_theta[i]  || theta_deg < min_theta[i]) return false; 
       var between = phi_deg <= max_phi[i] && phi_deg >= min_phi[i]; 
       return  (between && min_phi[i] < max_phi[i]) || (!between && min_phi[i] > max_phi[i]);
    }
  );


}


/** Computes a coherent sum of graphs with the given time delays */

RF.coherentSum = function( graphs, times ) 
{

  //set up the grid, use min / max and dt of first graph 
  if (!(graphs.length >0) || graphs.length!=times.length) return null; 
  var min = graphs[0].fX[0] - times[0]; 
  var max = graphs[0].fX[graphs[0].fNpoints-1] - times[0]; 
  var dt = (max-min)/(graphs[0].fNpoints -1); 

  for (var i = 1; i < N; i++) 
  {
    if (graphs[i].fX[0] -times[i]< min) min = graphs[i].fX[0]-times[i]; 
    if (graphs[i].fX[graphs[i].fNpoints-1] > max) max = graphs[i].fX[graphs[i].fNpoints-1] - times[i]; 
  }

  var N = (max-min)/dt+1; 

  var x = RF.range(min,N,dt); 
  var y= new Float32Array(N); 

  for (var i = 0; i < N; i++) 
  {
    for (var j = 0; j < graphs.length; j++)
    {
      y[i] += evalEven(graphs[j], x[i]-times[j]); 
    }
  }

  var ans = ROOT.CreateTGraph(N,x,y); 
  ans.fTitle = "Coherent Sum";
  ans.fXaxis.fTitle = "time"; 
  return ans; 

}


RF.Spectrogram = function(title, ntime, tmin, tmax, nfreq, fmin, fmax) 
{
  this.hist = JSROOT.CreateHistogram("TH2D", ntime, nfreq); 
  this.sum = JSROOT.CreateHistogram("TH2D", ntime, nfreq); 
  this.norm = JSROOT.CreateHistogram("TH1F", ntime); 
  this.tmin = tmin; 
  this.dt = (this.tmax-this.tmin)/(this.ntime-1); 
  this.nt = ntime; 
  this.nf = nfreq;
  this.df = (fmax-fmin)/(nfreq);
  this.fmin  = fmin; 

  this.norm.fXaxis.fXmin = tmin; 
  this.norm.fXaxis.fXmax = tmax; 
  this.norm.fXaxis.fTimeDisplay = 1; 

  this.hist.fTitle = title; 
  this.hist.fXaxis.fXmin = tmin; 
  this.hist.fXaxis.fXmax = tmax; 
  this.hist.fXaxis.fTimeDisplay = 1; 
  this.hist.fYaxis.fXmin = fmin; 
  this.hist.fYaxis.fXmax = fmax; 

  this.sum.fTitle = title; 
  this.sum.fXaxis.fXmin = tmin; 
  this.sum.fXaxis.fXmax = tmax; 
  this.sum.fXaxis.fTimeDisplay = 1; 
  this.sum.fYaxis.fXmin = fmin; 
  this.sum.fYaxis.fXmax = fmax; 






  this.addY = function(y,t) 
  {
    var Y = RF.doFFT(y); 
    var N = y.length; 
    var P = new Float32Array(N/2+1); 

    for (var i = 0; i <N/2+1; i++)
    {
      P[i] = (Y[2*i]*Y[2*i] + Y[2*i+1]*Y[2*i+1]) / N; 
      if (i > 0 || i < N/2) P[i] *=2; 
    }


    for (var i = 0; i < this.nf; i++) 
    {
      var f = this.fmin + (i+0.5) * this.df; 
      this.sum.Fill(t,f, P[i]); 
    }
    this.norm.Fill(t); 
  }


  this.finalize = function() 
  {

    for (var i = 0; i < this.nt; i++) 
    {
      for (var j = 0; j < this.nf; j++) 
      {
        var ibin = (this.nt+2) * (j+1) + i+1;
//        console.log(i,j,ibin, this.hist.getBinContent(ibin), this.norm.getBinContent(i+1)); 
//        this.hist.setBinContent(ibin, 10 * Math.log10(this.hist.getBinContent(ibin)/ this.norm.getBinContent(i+1))); 
        var n = this.norm.fArray[i+1]; 
        var s = this.sum.fArray[ibin]; 
        var r = s/n; 
        var dbish = n  == 0 ? -20 : 10*Math.log10(r); 
        if (dbish < -20) dbish = -20; 
        this.hist.setBinContent(ibin,  dbish); 
      }
    }
    this.hist.fXaxis.fTitle = "time"; 
    this.hist.fYaxis.fTitle = "freq"; 
    this.hist.fZaxis.fTitle = "dBish"; 
  }

}



/** This sets up an interferometric map (which will be a TH2) 
 *  with times defined by the mapper. The delta t's and antennas used
 *  are calculated once at the beginning, which is efficient when one keeps
 *  reusing the same map
 *
 * */ 
RF.InterferometricMap = function ( nx, xmin, xmax, ny, ymin,ymax, mapper) 
{

  /* make the histogram */ 
  this.hist = JSROOT.CreateHistogram("TH2F", nx,ny); 
  this.hist.fXaxis.fXmin = xmin; 
  this.hist.fXaxis.fXmax = xmax; 
  this.hist.fYaxis.fXmin = ymin; 
  this.hist.fYaxis.fXmax = ymax; 

  this.soln = RF.createArray(nx,ny); 

  this.usepair = RF.createArray(nx,ny);
  
  this.nx = nx; 
  this.ny = ny; 
  this.ymin = ymin; 
  this.ymax = ymax; 
  this.xmin = xmin; 
  this.xmax = xmax; 
  this.dx = (xmax-xmin)/nx; 
  this.dy = (ymax-ymin)/ny; 
  this.nant = mapper.nant;
  this.xcorrs = RF.createArray(this.nant, this.nant);

  this.cutoff = 0; 
  this.upsample = 4; 
  this.is_init = false; 

  this.init = function() 
  {
    if (this.is_init) return; 

    for (var iant = 0; iant < mapper.nant; iant++)
    {
      for (var jant = 0; jant < mapper.nant; jant++)
      {
        this.usepair[iant][jant] = mapper.usePair(iant,jant);
      }
    }


    for (var ix = 0; ix <this.nx; ix++) 
    {
      for (var iy = 0; iy < this.ny; iy++) 
      {
        var x = (ix+0.5)*this.dx + xmin; 
        var y = (iy+0.5)*this.dy + ymin; 

        this.soln[ix][iy] = []; 

        for (var iant  = 0; iant < mapper.nant; iant++)
        {
          if (!mapper.canUse(iant,x,y)) continue; 

          for (var jant = iant+1; jant < mapper.nant; jant++) 
          {
            if (!mapper.canUse(jant,x,y)) continue; 

            if (!this.usepair[iant][jant]) continue;  

            
            var this_soln = {}; 
            this_soln.i = iant; 
            this_soln.j = jant; 
            this_soln.dt = mapper.deltaTs(iant,jant,x,y);
            this.soln[ix][iy].push(this_soln); 
          }
        }
      }
    }
    this.is_init = true; 
  }

  this.setTitle = function (title, xtitle, ytitle) 
  {
    this.hist.fTitle = title;
    this.hist.fXaxis.fTitle = xtitle;
    this.hist.fYaxis.fTitle = ytitle;
  }


  this.deltaTHist = function (i,j)
  {
    this.init(); 
   
    if (i > j) 
    {
      var tmp  = i; 
      i = j; 
      j = tmp; 
    }

    var h = JSROOT.CreateHistogram("TH2F", nx,ny); 
    h.fXaxis.fXmin = xmin; 
    h.fXaxis.fXmax = xmax; 
    h.fYaxis.fXmin = ymin; 
    h.fYaxis.fXmax = ymax; 
    h.fXaxis.fTitle = "azimuth (degrees)"; 
    h.fYaxis.fTitle = "elevation (degrees)"; 
    h.fYaxis.fTitle = "elevation (degrees)"; 
    h.fTitle = "#Delta T (" + i + "," + j+")"; 

    for (var ix = 0; ix < this.nx; ix++)
    {
      for (var iy = 0; iy < this.ny; iy++)
      {
        var val = 0; 
        for (var ipair = 0; ipair < this.soln[ix][iy].length; ipair++)
        {
          if (this.soln[ix][iy][ipair].i == i && this.soln[ix][iy][ipair].j == j)
          {
            val = this.soln[ix][iy][ipair].dt; 
          }
        }
         var ibin = (this.nx+2) * (iy+1) + ix+1;
         h.setBinContent(ibin, val); 
      }
    }

    return h; 
  }

  this.compute = function(channels, reverse_sign = true) 
  {
    this.init(); 

    this.xcorrs = RF.createArray(this.nant, this.nant); 
    this.channels = channels; 

    for (var iant = 0; iant < this.nant; iant++) 
    {
      if (channels[iant] == null) continue;
      for (var jant = iant+1; jant < this.nant; jant++) 
      {
        if (channels[jant] == null) continue;
        if (this.usepair[iant][jant]) 
        {
          this.xcorrs[iant][jant] = RF.crossCorrelation(channels[iant], channels[jant],true,this.upsample,null,this.cutoff); 
        }
      }
    }

    for (var ix = 0; ix < this.nx; ix++) 
    {
      for (var iy = 0; iy < this.ny; iy++) 
      {
        var sum = 0; 
        var N = this.soln[ix][iy].length; 

        var norm = N; 
        for (var ipair = 0; ipair < N; ipair++) 
        {
          var soln = this.soln[ix][iy][ipair]; 
          if (channels[soln.i] == null || channels[soln.j] == null) 
          {
            norm--; 
            continue; 
          }
          sum += RF.evalEven(this.xcorrs[soln.i][soln.j], reverse_sign ? -soln.dt: soln.dt); 
        }

//        console.log(ix,iy, sum,norm);
        var ibin = (this.nx+2) * (iy+1) + ix+1;
        this.hist.setBinContent(ibin, norm? sum/norm : 0); 
      }
    }
  }

}



/** Perform an IIR filter (in-place) on a TGraph (g) 
 * using coefficients b and a
 * */
RF.IIRFilter = function(g, b,a) 
{

  if (a == null || a.length == 0) a = [1]; 
  if (b == null || b.length == 0) b = [1]; 

  var yNew = new Float32Array(g.fNpoints); 

  var inv = 1./a[0]; 
  for (var i = 0; i < g.fNpoints; i++) 
  {
    var val = 0; 

    for (var j = 0; j < Math.min(i+1,b.length); j++)
    {
      val += b[j] * g.fY[i-j]; 
    }

    for (var k = 1; k < Math.min(i+1,a.length); k++) 
    {
      val -= a[k] * yNew[i-k]; 
    }

    yNew[i] = val *inv; 

  }

  for (var i = 0; i < g.fNpoints;i++) 
  {
    g.fY[i] = yNew[i]; 
  }
}

