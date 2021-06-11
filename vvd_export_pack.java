import ij.*;
import ij.io.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import java.awt.image.*;
import ij.plugin.*;
import ij.plugin.filter.*;
import ij.plugin.frame.*;

import java.io.*;
import java.nio.*;
import java.util.*;
import java.util.zip.*;
import java.util.concurrent.atomic.AtomicInteger;

import javax.imageio.*;
import javax.imageio.stream.*;
import javax.xml.parsers.*;
import javax.xml.transform.*;
import javax.xml.transform.dom.*;
import javax.xml.transform.stream.*;

import org.apache.xml.serializer.OutputPropertiesFactory;
import org.w3c.dom.*;



public class vvd_export_pack implements PlugIn {
	String basename;
	String directory;
	ArrayList<String> lvImgTitle;
	int bw = 128, bh = 128, bd = 64, lv = 1;
	ArrayList<Integer> bwlist = new ArrayList<Integer>();
	ArrayList<Integer> bhlist = new ArrayList<Integer>();
	ArrayList<Integer> bdlist = new ArrayList<Integer>();
	int thread_num_ = (int)Prefs.get("thread_num.int",4);
	String filetype = (String)Prefs.get("filetype.string", "RAW");
	int jpeg_quality = (int)Prefs.get("jpeg_quality.int",FileSaver.DEFAULT_JPEG_QUALITY);


	class Brick {
		public int x_, y_, z_;
		public int w_, h_, d_;
		public long offset_, size_;
		public double tx0_, ty0_, tz0_, tx1_, ty1_, tz1_;
		public double bx0_, by0_, bz0_, bx1_, by1_, bz1_;
		Brick(int x, int y, int z, int w, int h, int d, long offset, long size,
			  double tx0, double ty0, double tz0, double tx1, double ty1, double tz1,
			  double bx0, double by0, double bz0, double bx1, double by1, double bz1){
			x_ = x; y_ = y; z_ = z;
			w_ = w; h_ = h; d_ = d;
			offset_ = offset; size_ = size;
			tx0_ = tx0; ty0_ = ty0; tz0_ = tz0; tx1_ = tx1; ty1_ = ty1; tz1_ = tz1;
			bx0_ = bx0; by0_ = by0; bz0_ = bz0; bx1_ = bx1; by1_ = by1; bz1_ = bz1;
		}
	}

	public void run(String arg) {
		if (IJ.versionLessThan("1.49d")) return;
		
		if (!showDialog()) return;
		if (!showDialog2()) return;
		
		SaveDialog sd = new SaveDialog("Save as Bricks...", "", "");
		basename = sd.getFileName();
		directory = sd.getDirectory();
		if(basename == null || directory == null) return;
		
		build_bricks();
	}
	
	private boolean showDialog() {
		String[] types = {"RAW", "JPEG", "ZLIB"};
		GenericDialog gd = new GenericDialog("Generate Bricks");
		gd.addNumericField("BrickWidth  (default)",  bw, 0);
		gd.addNumericField("BrickHeight (default)", bh, 0);
		gd.addNumericField("BrickDepth  (default)",  bd, 0);
		gd.addNumericField("Levels", lv, 0);
		gd.addChoice("FileType", types, types[0]);
		gd.addNumericField("Thread", thread_num_, 0);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		bw = (int)gd.getNextNumber();
		bh = (int)gd.getNextNumber();
		bd = (int)gd.getNextNumber();
		lv = (int)gd.getNextNumber();
		filetype = types[gd.getNextChoiceIndex()];
		thread_num_ = (int)gd.getNextNumber();
		if (thread_num_ <= 0) thread_num_ = 1;
		
		Prefs.set("thread_num.int",thread_num_);
		return true;
	}

	private boolean showDialog2() {
		
		int[] wlist = WindowManager.getIDList();
		if(wlist == null)return false;

		String[] titles = new String[wlist.length+1];
		for(int i = 0; i < wlist.length; i++) titles[i] = "";
		
		int tnum = 0;
		for(int i = 0; i < wlist.length; i++){
			ImagePlus imp = WindowManager.getImage(wlist[i]);
			if(imp != null){
				titles[tnum] = imp.getTitle();
				tnum++;
			}
		}
		titles[tnum] = "none";
		
		GenericDialog gd = new GenericDialog("Files");
		if(filetype == "JPEG")gd.addNumericField("JPEG quality", jpeg_quality, 0);
		if(lv <= 0) return false;
		Font font = new Font("Dialog", Font.BOLD, 13);
		for(int i = 0; i < lv; i++) { 
			if(i == 0) gd.setInsets(0,0,0);
			if(i > 0) gd.setInsets(5,0,0);
			gd.addMessage("Level:" + String.valueOf(i), font);
			gd.addChoice("Image"+ String.valueOf(i)+":", titles, i < tnum ? titles[i] : titles[tnum] );
			gd.addStringField("Bricksize"+ String.valueOf(i)+" (w h d)",  bw + " " + bh + " " + bd);
		}
		gd.showDialog();
		if (gd.wasCanceled()) return false;

		if(filetype == "JPEG")jpeg_quality = (int)gd.getNextNumber();
		if(jpeg_quality > 100)jpeg_quality = 100;
		if(jpeg_quality < 0)  jpeg_quality = 0;
		
		lvImgTitle = new ArrayList<String>();
		for(int i = 0; i < lv; i++){
			int id = gd.getNextChoiceIndex();
			String bsizetxt = gd.getNextString();
			IJ.log(bsizetxt);
			if(bsizetxt != null && id < tnum){
				String[] bwhd = bsizetxt.split(" ");
				if(bwhd.length == 3){
					bwlist.add(Integer.parseInt(bwhd[0]));
					bhlist.add(Integer.parseInt(bwhd[1]));
					bdlist.add(Integer.parseInt(bwhd[2]));
				}
				else{
					bwlist.add(bw);
					bhlist.add(bh);
					bdlist.add(bd);
				}
				lvImgTitle.add(titles[id]);
			}
		}
		for(int i = 0; i < lvImgTitle.size(); i++) IJ.log(lvImgTitle.get(i));

		if(lvImgTitle.isEmpty()) return false;
		
		return true;
	}

	public void build_bricks(){

		ImagePlus imp;
		ImagePlus orgimp;
		ImageStack stack;
		FileInfo finfo;
				
		if(lvImgTitle.isEmpty()) return;
		orgimp = WindowManager.getImage(lvImgTitle.get(0));
		imp = orgimp;
		
		finfo = imp.getFileInfo();
		if(finfo == null) return;

		int[] dims = imp.getDimensions();
		int imageW = dims[0];
		int imageH = dims[1];
		int nCh    = dims[2];
		int imageD = dims[3];
		int nFrame = dims[4];
		int bdepth = imp.getBitDepth();
		double xspc = finfo.pixelWidth;
		double yspc = finfo.pixelHeight;
		double zspc = finfo.pixelDepth;

		int orgW = imageW;
		int orgH = imageH;
		int orgD = imageD;
		double orgxspc = xspc;
		double orgyspc = yspc;
		double orgzspc = zspc;

		lv = lvImgTitle.size();
		if(filetype == "JPEG"){
			for(int l = 0; l < lv; l++){
				if (WindowManager.getImage(lvImgTitle.get(l)).getBitDepth() != 8) {
					IJ.error("ALL IMAGES MUST BE 8BIT GLAYSCALE");
					return;
				}
			}
		}

		Document doc = newXMLDocument();
		Element root = doc.createElement("BRK");
		root.setAttribute("version", "1.0");
		root.setAttribute("nLevel", String.valueOf(lv));
		root.setAttribute("nChannel", String.valueOf(nCh));
		root.setAttribute("nFrame", String.valueOf(nFrame));
		doc.appendChild(root);

		for(int l = 0; l < lv; l++){
			
			int[] dims2 = imp.getDimensions();
			IJ.log("W: " + String.valueOf(dims2[0]) +
				  " H: " + String.valueOf(dims2[1]) +
				  " C: " + String.valueOf(dims2[2]) +
				  " D: " + String.valueOf(dims2[3]) +
				  " T: " + String.valueOf(dims2[4]) +
				  " b: " + String.valueOf(bdepth)   );

			bw = bwlist.get(l).intValue();
			bh = bhlist.get(l).intValue();
			bd = bdlist.get(l).intValue();

			if (l == lv - 1) {
				bw = imageW;
				bh = imageH;
				bd = imageD;
			}

			boolean force_pow2 = false;
			if(IsPowerOf2(bw) && IsPowerOf2(bh) && IsPowerOf2(bd)) force_pow2 = true;
			
			if(force_pow2){
				//force pow2
				if(Pow2(bw) > bw) bw = Pow2(bw)/2;
				if(Pow2(bh) > bh) bh = Pow2(bh)/2;
				if(Pow2(bd) > bd) bd = Pow2(bd)/2;
			}

			if(bw > imageW) bw = (Pow2(imageW) == imageW) ? imageW : Pow2(imageW)/2;
			if(bh > imageH) bh = (Pow2(imageH) == imageH) ? imageH : Pow2(imageH)/2;
			if(bd > imageD) bd = (Pow2(imageD) == imageD) ? imageD : Pow2(imageD)/2;

			if(bw <= 1 || bh <= 1 || bd <= 1) break;

			if(filetype == "JPEG" && (bw < 8 || bh < 8)) break;
		
			Element lvnode = doc.createElement("Level");
			lvnode.setAttribute("lv", String.valueOf(l));
			lvnode.setAttribute("imageW", String.valueOf(imageW));
			lvnode.setAttribute("imageH", String.valueOf(imageH));
			lvnode.setAttribute("imageD", String.valueOf(imageD));
			lvnode.setAttribute("xspc", String.valueOf(xspc));
			lvnode.setAttribute("yspc", String.valueOf(yspc));
			lvnode.setAttribute("zspc", String.valueOf(zspc));
			lvnode.setAttribute("bitDepth", String.valueOf(bdepth));
			lvnode.setAttribute("FileType", filetype);
			root.appendChild(lvnode);

			Element brksnode = doc.createElement("Bricks");
			brksnode.setAttribute("brick_baseW", String.valueOf(bw));
			brksnode.setAttribute("brick_baseH", String.valueOf(bh));
			brksnode.setAttribute("brick_baseD", String.valueOf(bd));
			lvnode.appendChild(brksnode);

			final ArrayList<Brick> bricks = new ArrayList<Brick>();
			int mw, mh, md, mw2, mh2, md2;
			double tx0, ty0, tz0, tx1, ty1, tz1;
			double bx0, by0, bz0, bx1, by1, bz1;
			ArrayList<Integer> xy_brk_count = new ArrayList<Integer>();
			for(int k = 0; k < imageD; k += bd){
				int xy_count = 0;
				if(k > 0) k--;
				for(int j = 0; j < imageH; j += bh){
					if(j > 0) j--;
					for(int i = 0; i < imageW; i += bw){
						if(i > 0) i--;
						mw = Math.min(bw, imageW - i);
						mh = Math.min(bh, imageH - j);
						md = Math.min(bd, imageD - k);

						if(force_pow2){
							mw2 = Pow2(mw);
							mh2 = Pow2(mh);
							md2 = Pow2(md);
						}
						else{
							mw2 = mw;
							mh2 = mh;
							md2 = md;
						}

						if (filetype == "JPEG") {
							if (mw2 < 8) mw2 = 8;
							if (mh2 < 8) mh2 = 8;
						}

						tx0 = i == 0 ? 0.0d : ((mw2 - mw + 0.5d) / mw2);
						ty0 = j == 0 ? 0.0d : ((mh2 - mh + 0.5d) / mh2);
						tz0 = k == 0 ? 0.0d : ((md2 - md + 0.5d) / md2);

						tx1 = 1.0d - 0.5d / mw2;
						if (mw < bw) tx1 = 1.0d;
						if (imageW - i == bw) tx1 = 1.0d;
	
						ty1 = 1.0d - 0.5d / mh2;
						if (mh < bh) ty1 = 1.0d;
						if (imageH - j == bh) ty1 = 1.0d;

						tz1 = 1.0d - 0.5d / md2;
						if (md < bd) tz1 = 1.0d;
						if (imageD - k == bd) tz1 = 1.0d;

						bx0 = i == 0 ? 0.0d : (i + 0.5d)/(double)imageW;
						by0 = j == 0 ? 0.0d : (j + 0.5d)/(double)imageH;
						bz0 = k == 0 ? 0.0d : (k + 0.5d)/(double)imageD;

						bx1 = Math.min((i + bw - 0.5d)/(double)imageW, 1.0d);
						if (imageW - i == bw) bx1 = 1.0d;
	
						by1 = Math.min((j + bh - 0.5d)/(double)imageH, 1.0d);
						if (imageH - j == bh) by1 = 1.0d;

						bz1 = Math.min((k + bd - 0.5d)/(double)imageD, 1.0d);
						if (imageD - k == bd) bz1 = 1.0d;

						int x, y, z;
						x = i - (mw2 - mw);
						y = j - (mh2 - mh);
						z = k - (md2 - md);
						bricks.add(new Brick(x, y, z, mw2, mh2, md2, 0, 0, tx0, ty0, tz0, tx1, ty1, tz1, bx0, by0, bz0, bx1, by1, bz1));
						xy_count++;
					}
				}
				xy_brk_count.add(xy_count);
			}

			Element fsnode = doc.createElement("Files");
			lvnode.appendChild(fsnode);

			stack = imp.getStack();

			int totalbricknum = nFrame * nCh * bricks.size();
			int curbricknum = 0;
			for(int f = 0; f < nFrame; f++){
				for(int ch = 0; ch < nCh; ch++){
					int first_bid = 0;
					int st_z = 0;
					int ed_z = 0;
					LinkedList<ImageProcessor> iplist = new LinkedList<ImageProcessor>();
					
					long bytecount = 0;
					long filecount = 0;
					String base_dataname = basename + "_Lv" + String.valueOf(l) +
													  "_Ch" + String.valueOf(ch)+
													  "_Fr" + String.valueOf(f);
					String current_dataname = base_dataname+"_data"+filecount;
					
					for (int ll = 0; ll < xy_brk_count.size(); ll++) {
						Brick b_first = bricks.get(first_bid);

						if (ll == 0) {
							st_z = b_first.z_;
							ed_z = b_first.z_ + b_first.d_;
							for(int s = st_z; s < ed_z; s++)
								iplist.add( stack.getProcessor(imp.getStackIndex(ch+1, s+1, f+1)) );
						}
						else if(ed_z > b_first.z_ || st_z < b_first.z_+b_first.d_){
							if(b_first.z_ > st_z){
								for(int s = 0; s < b_first.z_-st_z; s++)iplist.pollFirst();
								st_z = b_first.z_;
							}
							else if(b_first.z_ < st_z){
								IJ.log("warning");
								for(int s = st_z-1; s > b_first.z_; s--)
									iplist.addFirst( stack.getProcessor(imp.getStackIndex(ch+1, s+1, f+1)) );
								st_z = b_first.z_;
							}
							
							if(b_first.z_+b_first.d_ > ed_z){
								for(int s = ed_z; s < b_first.z_+b_first.d_; s++)
									iplist.add( stack.getProcessor(imp.getStackIndex(ch+1, s+1, f+1)) );
									ed_z = b_first.z_+b_first.d_;
							}
							else if(b_first.z_+b_first.d_ < ed_z){
								IJ.log("warning");
								for(int s = 0; s < ed_z-(b_first.z_+b_first.d_); s++)iplist.pollLast();
								ed_z = b_first.z_+b_first.d_;
							}
						}
						else {
							IJ.log("warning");
							iplist.clear();
							st_z = b_first.z_;
							ed_z = b_first.z_+b_first.d_;
							for(int s = st_z; s < ed_z; s++)
							iplist.add( stack.getProcessor(imp.getStackIndex(ch+1, s+1, f+1)) );
						}

						if(iplist.size() != b_first.d_){
							IJ.log("Stack Error "+iplist.size()+" "+b_first.d_);
							return;
						}

						final int st_bid = first_bid;
						final int ed_bid = first_bid + xy_brk_count.get(ll) - 1;

						final AtomicInteger ai = new AtomicInteger(st_bid);
						final Thread[] threads = newThreadArray();

						Iterator<ImageProcessor> ipite = iplist.iterator();
						final ArrayList<ImageProcessor> iplist_final = new ArrayList<ImageProcessor>();
						while (ipite.hasNext()){
							iplist_final.add(ipite.next());
						}

						final int pitchY = imageW;
						final String fftype = filetype;
						final String fbasename = basename;
						final ArrayList<Element> elems = new ArrayList<Element>();
						for (int e = 0; e < xy_brk_count.get(ll); e++)
							elems.add(doc.createElement("File"));

						final byte[][] data_array = new byte[xy_brk_count.get(ll)][]; 

						final int fbd = bdepth;
						final int flv = l;
						final int fch = ch;
						final int ffr = f;

						final String fdataname = current_dataname;
		
						for (int ithread = 0; ithread < threads.length; ithread++) {
							// Concurrently run in as many threads as CPUs
							threads[ithread] = new Thread() {
								
								{ setPriority(Thread.NORM_PRIORITY); }
				
								public void run() {
					
									for (int i = ai.getAndIncrement(); i <= ed_bid; i = ai.getAndIncrement()) {
										Brick b = bricks.get(i);
										int bsize = 0;
										byte[] bdata = new byte[b.w_*b.h_*b.d_*fbd/8];

										for (int p = 0; p < iplist_final.size(); p++){
											ImageProcessor ip = iplist_final.get(p);
											ip.setRoi(b.x_, b.y_, b.w_, b.h_);
											if(fbd == 8){
												byte[] data = (byte[])ip.getPixels();
												for (int yy = b.y_; yy < b.y_+b.h_; yy++) {
													int srcpos = yy*pitchY + b.x_;
													System.arraycopy(data, srcpos, bdata, bsize, b.w_);
													bsize += b.w_;
												}
												data = null;
											} else if (fbd == 16) {
												ByteBuffer buffer = ByteBuffer.allocate(b.w_*b.h_*fbd/8);
												buffer.order(ByteOrder.LITTLE_ENDIAN);
												short[] data = (short[])ip.getPixels();
												for (int yy = b.y_; yy < b.y_+b.h_; yy++) {
													for (int xx = b.x_; xx < b.x_+b.w_; xx++) {
														buffer.putShort((short)(data[yy*pitchY+xx]<<4));
													}
												}
												System.arraycopy(buffer.array(), 0, bdata, bsize, buffer.array().length);
												bsize += buffer.array().length;
												data = null;
												buffer = null;
											} else if (fbd == 32) {
												ByteBuffer buffer = ByteBuffer.allocate(b.w_*b.h_*fbd/8);
												buffer.order(ByteOrder.LITTLE_ENDIAN);
												float[] data = (float[])ip.getPixels();
												for (int yy = b.y_; yy < b.y_+b.h_; yy++) {
													for (int xx = b.x_; xx < b.x_+b.w_; xx++) {
														buffer.putFloat(data[yy*pitchY+xx]);
													}
												}
												System.arraycopy(buffer.array(), 0, bdata, bsize, buffer.array().length);
												bsize += buffer.array().length;
												data = null;
												buffer = null;
											}
										}

										int datasize = bdata.length;

										if (fftype.equals("RAW")) {
											int dummy = -1;
											//do nothing
										}
										if (fftype.equals("JPEG") && fbd == 8) {
											try {
												DataBufferByte db = new DataBufferByte(bdata, datasize);
												Raster raster = Raster.createPackedRaster(db, b.w_, b.h_*b.d_, 8, null);
												BufferedImage img = new BufferedImage(b.w_, b.h_*b.d_, BufferedImage.TYPE_BYTE_GRAY);
												img.setData(raster);
												ByteArrayOutputStream baos = new ByteArrayOutputStream();
												ImageOutputStream ios = ImageIO.createImageOutputStream(baos);
												String format = "jpg";
												Iterator<javax.imageio.ImageWriter> iter = ImageIO.getImageWritersByFormatName("jpeg");
												javax.imageio.ImageWriter writer = iter.next();
												ImageWriteParam iwp = writer.getDefaultWriteParam();
												iwp.setCompressionMode(ImageWriteParam.MODE_EXPLICIT);
												iwp.setCompressionQuality((float)jpeg_quality*0.01f);
												writer.setOutput(ios);
												writer.write(null, new IIOImage(img,null,null), iwp);
												bdata = baos.toByteArray();
												datasize = bdata.length;
											} catch (IOException e) {
												e.printStackTrace();
												return;
											}
										}
										if (fftype.equals("ZLIB")) {
											byte[] tmpdata = new byte[b.w_*b.h_*b.d_*fbd/8];
											Deflater compresser = new Deflater(Deflater.BEST_SPEED);
											compresser.setInput(bdata);
											//compresser.setLevel(Deflater.DEFAULT_COMPRESSION);
											//compresser.setStrategy(Deflater.DEFAULT_STRATEGY);
											compresser.finish();
											datasize = compresser.deflate(tmpdata);
											bdata = tmpdata;
											compresser.end();
										}

										if (bdata.length != datasize)
										{
											byte[] tmpdata2 = new byte[datasize];
											System.arraycopy(bdata, 0, tmpdata2, 0, datasize);
											bdata = tmpdata2;
										}

										data_array[i-st_bid] = bdata;

										Element filenode = elems.get(i - st_bid);
										filenode.setAttribute("filename", fdataname);
										filenode.setAttribute("channel", String.valueOf(fch));
										filenode.setAttribute("frame", String.valueOf(ffr));
										filenode.setAttribute("brickID", String.valueOf(i));
										filenode.setAttribute("filetype", String.valueOf(fftype));

									}//	for (int i = ai.getAndIncrement(); i < names.length;
							}};//threads[ithread] = new Thread() {
						}//	for (int ithread = 0; ithread < threads.length; ithread++) {
						
						startAndJoin(threads);
						
						BufferedOutputStream fis = null;
						try {
							File file = new File(directory + current_dataname);
							fis = new BufferedOutputStream(new FileOutputStream(file, bytecount == 0 ? false : true));
							for (int d = 0; d < data_array.length; d++) {
								long offset = bytecount;
								long datasize = data_array[d].length;
								elems.get(d).setAttribute("offset", String.valueOf(offset));
								elems.get(d).setAttribute("datasize", String.valueOf(datasize));
								fsnode.appendChild(elems.get(d));
								fis.write(data_array[d]);
								bytecount += datasize;
							}
						} catch (IOException e) {
							e.printStackTrace();
							return;
						} finally {
							try {
								if (fis != null) fis.close();
							} catch (IOException e) {
								e.printStackTrace();
								return;
							}
						}
						
						first_bid += xy_brk_count.get(ll);
						System.gc();
					}
					iplist.clear();
					System.gc();

				}//for(int ch = 0; ch < nCh; ch++)
			}//for(int f = 0; f < nFrame; f++)

			for(int i = 0; i < bricks.size(); i++){
				Brick b = bricks.get(i);
				Element bricknode = doc.createElement("Brick");
				bricknode.setAttribute("id", String.valueOf(i));
				bricknode.setAttribute("st_x", String.valueOf(b.x_));
				bricknode.setAttribute("st_y", String.valueOf(b.y_));
				bricknode.setAttribute("st_z", String.valueOf(b.z_));
				bricknode.setAttribute("width", String.valueOf(b.w_));
				bricknode.setAttribute("height", String.valueOf(b.h_));
				bricknode.setAttribute("depth", String.valueOf(b.d_));
				brksnode.appendChild(bricknode);

				Element tboxnode = doc.createElement("tbox");
				tboxnode.setAttribute("x0", String.valueOf(b.tx0_));
				tboxnode.setAttribute("y0", String.valueOf(b.ty0_));
				tboxnode.setAttribute("z0", String.valueOf(b.tz0_));
				tboxnode.setAttribute("x1", String.valueOf(b.tx1_));
				tboxnode.setAttribute("y1", String.valueOf(b.ty1_));
				tboxnode.setAttribute("z1", String.valueOf(b.tz1_));
				bricknode.appendChild(tboxnode);

				Element bboxnode = doc.createElement("bbox");
				bboxnode.setAttribute("x0", String.valueOf(b.bx0_));
				bboxnode.setAttribute("y0", String.valueOf(b.by0_));
				bboxnode.setAttribute("z0", String.valueOf(b.bz0_));
				bboxnode.setAttribute("x1", String.valueOf(b.bx1_));
				bboxnode.setAttribute("y1", String.valueOf(b.by1_));
				bboxnode.setAttribute("z1", String.valueOf(b.bz1_));
				bricknode.appendChild(bboxnode);
			}
						

			if(l < lv-1){
				imp = WindowManager.getImage(lvImgTitle.get(l+1));
				int[] newdims = imp.getDimensions();
				imageW = newdims[0];
				imageH = newdims[1];
				imageD = newdims[3];
				xspc = orgxspc * ((double)orgW / (double)imageW);
				yspc = orgyspc * ((double)orgH / (double)imageH);
				zspc = orgzspc * ((double)orgD / (double)imageD);
				bdepth = imp.getBitDepth();
			}
			
		}

		try{
			File newXMLfile = new File(directory + basename + ".vvd");
			newXMLfile.createNewFile();
			writeXML(newXMLfile, doc);
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}

	
	}


	public static Document newXMLDocument() {
		DocumentBuilder dbuilder = null;
		try {
			dbuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
		} catch (ParserConfigurationException e) {
			e.printStackTrace();
			return null;
		}
		return dbuilder.newDocument();
	}

	public static boolean writeXML(File file, Document document) {

		javax.xml.transform.Transformer transformer = null;
		try {
			transformer = TransformerFactory.newInstance().newTransformer();
		} catch (TransformerConfigurationException e) {
			e.printStackTrace();
			return false;
		}

		transformer.setOutputProperty("indent", "yes");
		transformer.setOutputProperty(OutputPropertiesFactory.S_KEY_INDENT_AMOUNT, "2");
		transformer.setOutputProperty("encoding", "UTF-8");

		try {
			transformer.transform(new DOMSource(document), new StreamResult(file));
		} catch (TransformerException e) {
			e.printStackTrace();
			return false;
		}

		return true;
	}

	// Fast way to check for power of two
	public static boolean IsPowerOf2(int n)
	{
		return (n & (n-1)) == 0;
	}

	// Returns a number Greater Than or Equal to dim
	// that is an exact power of 2
	// Used for determining what size of texture to
	// allocate to store an image
	public static int Pow2(int dim)
	{
		if (IsPowerOf2(dim)) return dim;
		int val = 4;
		while (val < dim) val = val << 1;
		return val;
	}

	/** Create a Thread[] array as large as the number of processors available.
		* From Stephan Preibisch's Multithreading.java class. See:
		* http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
		*/
	private Thread[] newThreadArray() {
		int n_cpus = Runtime.getRuntime().availableProcessors();
		if (n_cpus > thread_num_) n_cpus = thread_num_;
		if (n_cpus <= 0) n_cpus = 1;
		return new Thread[n_cpus];
	}
	
	/** Start all given threads and wait on each of them until all are done.
		* From Stephan Preibisch's Multithreading.java class. See:
		* http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
		*/
	public static void startAndJoin(Thread[] threads)
	{
		for (int ithread = 0; ithread < threads.length; ++ithread)
		{
			threads[ithread].setPriority(Thread.NORM_PRIORITY);
			threads[ithread].start();
		}
		
		try
		{   
			for (int ithread = 0; ithread < threads.length; ++ithread)
			threads[ithread].join();
		} catch (InterruptedException ie)
		{
			throw new RuntimeException(ie);
		}
	}

}
