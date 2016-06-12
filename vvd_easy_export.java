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

import javax.imageio.*;
import javax.imageio.stream.*;
import javax.xml.parsers.*;
import javax.xml.transform.*;
import javax.xml.transform.dom.*;
import javax.xml.transform.stream.*;

import org.apache.xml.serializer.OutputPropertiesFactory;
import org.w3c.dom.*;



public class vvd_easy_export implements PlugIn {
	String basename;
	String directory;
	ArrayList<String> lvImgTitle;
	int bw = 128, bh = 128, bd = 64, lv = 1;
	String filetype = "RAW";
	int jpeg_quality = FileSaver.DEFAULT_JPEG_QUALITY;
	ArrayList<Integer> bwlist = new ArrayList<Integer>();
	ArrayList<Integer> bhlist = new ArrayList<Integer>();
	ArrayList<Integer> bdlist = new ArrayList<Integer>();

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
				
		SaveDialog sd = new SaveDialog("Save as Bricks...", "", "");
		basename = sd.getFileName();
		directory = sd.getDirectory();
		if(basename == null || directory == null) return;
		
		build_bricks();
	}
	
	private boolean showDialog() {
		String[] types = {"RAW", "JPEG", "ZLIB"};
		GenericDialog gd = new GenericDialog("Generate Bricks");
		gd.addChoice("FileType", types, types[0]);
		gd.addNumericField("JPEG quality", jpeg_quality, 0);


		int[] wlist = WindowManager.getIDList();
		if(wlist == null)return false;

		String[] titles = new String[wlist.length];
		for(int i = 0; i < wlist.length; i++) titles[i] = "";
		
		int tnum = 0;
		for(int i = 0; i < wlist.length; i++){
			ImagePlus imp = WindowManager.getImage(wlist[i]);
			if(imp != null){
				titles[tnum] = imp.getTitle();
				tnum++;
			}
		}
		gd.addChoice("Source image: ", titles, titles[0]);
		
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		
		filetype = types[gd.getNextChoiceIndex()];
		if(filetype == "JPEG")jpeg_quality = (int)gd.getNextNumber();
		if(jpeg_quality > 100)jpeg_quality = 100;
		if(jpeg_quality < 0)  jpeg_quality = 0;

		int id = gd.getNextChoiceIndex();
		lvImgTitle = new ArrayList<String>();
		lvImgTitle.add(titles[id]);
		
		return true;
	}

	public static int log2(int n){
    	if(n <= 0) return 0;
    	return 31 - Integer.numberOfLeadingZeros(n);
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
		double z_aspect = Math.max(xspc, yspc)/zspc;

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
					IJ.error("A SOURCE IMAGE MUST BE 8BIT GLAYSCALE");
					return;
				}
			}
		}

		//calculate levels
/*		int baseXY = 256;
		int baseZ = 256;

		if (z_aspect < 0.5) baseZ = 128;
		if (z_aspect > 2.0) baseXY = 128;
		if (z_aspect >= 0.5 && z_aspect < 1.0) baseZ = (int)(baseZ*z_aspect); 
		if (z_aspect > 1.0 && z_aspect <= 2.0) baseXY = (int)(baseXY/z_aspect); 

		IJ.log("Z_aspect: " + z_aspect);
		IJ.log("BaseXY: " + baseXY);
		IJ.log("BaseZ: " + baseZ);
*/

		int baseXY = 256;
		int baseZ = 128;
		int dbXY = Math.max(orgW, orgH)/baseXY;
		if (Math.max(orgW, orgH) % baseXY > 0) dbXY *= 2;
		int dbZ = orgD/baseZ;
		if (orgD % baseZ > 0) dbZ *= 2;
		lv = Math.max(log2(dbXY), log2(dbZ)) + 1;

		int ww = orgW;
		int hh = orgH;
		int dd = orgD;
		for (int l = 0; l < lv; l++){
			int bwnum = ww / baseXY; if (ww % baseXY > 0) bwnum++;
			int bhnum = hh / baseXY; if (hh % baseXY > 0) bhnum++;
			int bdnum = dd / baseZ;  if (dd % baseZ  > 0) bdnum++;

			if (bwnum % 2 == 0) bwnum++;
			if (bhnum % 2 == 0) bhnum++;
			if (bdnum % 2 == 0) bdnum++;

			int bw = (bwnum <= 1) ? ww : ww/bwnum+1+(ww%bwnum>0 ? 1 : 0);
			int bh = (bhnum <= 1) ? hh : hh/bhnum+1+(hh%bhnum>0 ? 1 : 0);
			int bd = (bdnum <= 1) ? dd : dd/bdnum+1+(dd%bdnum>0 ? 1 : 0);

			bwlist.add(bw);
			bhlist.add(bh);
			bdlist.add(bd);

			IJ.log("LEVEL: " + l);
			IJ.log("  width: " + ww);
			IJ.log("  hight: " + hh);
			IJ.log("  depth: " + dd);
			IJ.log("  bw: " + bw);
			IJ.log("  bh: " + bh);
			IJ.log("  bd: " + bd);
			
			int xyl2 = Math.max(ww, hh)/baseXY;
			if (Math.max(ww, hh) % baseXY > 0) xyl2 *= 2;
			if (lv-1-log2(xyl2) <= l) { ww /= 2; hh /= 2; }
			IJ.log("  xyl2: " + (lv-1-log2(xyl2)));
			
			int zl2 = dd/baseZ;
			if (dd % baseZ > 0) zl2 *= 2;
			if (lv-1-log2(zl2) <= l) dd /= 2;
			IJ.log("  zl2: " + (lv-1-log2(zl2)));

			if (l < lv-1) {
				lvImgTitle.add(lvImgTitle.get(0)+"_level"+(l+1));
				IJ.selectWindow(lvImgTitle.get(0));
				IJ.run("Scale...", "x=- y=- z=- width="+ww+" height="+hh+" depth="+dd+" interpolation=Bicubic average process create title="+lvImgTitle.get(l+1));
			}
		}
		
		for (int l = 0; l < lv; l++){
			IJ.log(lvImgTitle.get(l));
		}

		Document doc = newXMLDocument();
		Element root = doc.createElement("BRK");
		root.setAttribute("version", "1.0");
		root.setAttribute("nLevel", String.valueOf(lv));
		root.setAttribute("nChannel", String.valueOf(nCh));
		root.setAttribute("nFrame", String.valueOf(nFrame));
		doc.appendChild(root);

		for(int l = 0; l < lv; l++){
			IJ.showProgress(0.0);
			
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

			ArrayList<Brick> bricks = new ArrayList<Brick>();
			int mw, mh, md, mw2, mh2, md2;
			double tx0, ty0, tz0, tx1, ty1, tz1;
			double bx0, by0, bz0, bx1, by1, bz1;
			for(int k = 0; k < imageD; k += bd){
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
					}
				}
			}

			Element fsnode = doc.createElement("Files");
			lvnode.appendChild(fsnode);

			stack = imp.getStack();

			int totalbricknum = nFrame * nCh * bricks.size();
			int curbricknum = 0;
			for(int f = 0; f < nFrame; f++){
				for(int ch = 0; ch < nCh; ch++){
					int sizelimit = Math.max(50000000, bw*bh*bd*bdepth/8);
					int bytecount = 0;
					int filecount = 0;
					byte[] packed_data = new byte[sizelimit];
					String base_dataname = basename + "_Lv" + String.valueOf(l) +
													  "_Ch" + String.valueOf(ch)+
													  "_Fr" + String.valueOf(f);
					String current_dataname = base_dataname+"_data"+filecount;
					
					Brick b_first = bricks.get(0);
					if(b_first.z_ != 0) IJ.log("warning");
					int st_z = b_first.z_;
					int ed_z = b_first.z_ + b_first.d_;
					LinkedList<ImageProcessor> iplist = new LinkedList<ImageProcessor>();
					for(int s = st_z; s < ed_z; s++)
						iplist.add( stack.getProcessor(imp.getStackIndex(ch+1, s+1, f+1)) );

//					ImagePlus test;
//					ImageStack tsst;
//					test = NewImage.createByteImage("test", imageW, imageH, imageD, NewImage.FILL_BLACK);
//					tsst = test.getStack();
					for(int i = 0; i < bricks.size(); i++){
						Brick b = bricks.get(i);

						if(ed_z > b.z_ || st_z < b.z_+b.d_){
							if(b.z_ > st_z){
								for(int s = 0; s < b.z_-st_z; s++)iplist.pollFirst();
								st_z = b.z_;
							}
							else if(b.z_ < st_z){
								IJ.log("warning");
								for(int s = st_z-1; s > b.z_; s--)
									iplist.addFirst( stack.getProcessor(imp.getStackIndex(ch+1, s+1, f+1)) );
								st_z = b.z_;
							}
							
							if(b.z_+b.d_ > ed_z){
								for(int s = ed_z; s < b.z_+b.d_; s++)
									iplist.add( stack.getProcessor(imp.getStackIndex(ch+1, s+1, f+1)) );
									ed_z = b.z_+b.d_;
							}
							else if(b.z_+b.d_ < ed_z){
								IJ.log("warning");
								for(int s = 0; s < ed_z-(b.z_+b.d_); s++)iplist.pollLast();
								ed_z = b.z_+b.d_;
							}
						}
						else {
							IJ.log("warning");
							iplist.clear();
							st_z = b.z_;
							ed_z = b.z_+b.d_;
							for(int s = st_z; s < ed_z; s++)
							iplist.add( stack.getProcessor(imp.getStackIndex(ch+1, s+1, f+1)) );
						}

						if(iplist.size() != b.d_){
							IJ.log("Stack Error");
							return;
						}

//						int zz = st_z;

						int bsize = 0;
						byte[] bdata = new byte[b.w_*b.h_*b.d_*bdepth/8];
						Iterator<ImageProcessor> ipite = iplist.iterator();
						while (ipite.hasNext()){

//							ImageProcessor tsip = tsst.getProcessor(zz+1);
							
							ImageProcessor ip = ipite.next();
							ip.setRoi(b.x_, b.y_, b.w_, b.h_);
							if(bdepth == 8){
								byte[] data = (byte[])ip.crop().getPixels();
								System.arraycopy(data, 0, bdata, bsize, data.length);
								bsize += data.length;
							}
							else if(bdepth == 16){
								ByteBuffer buffer = ByteBuffer.allocate(b.w_*b.h_*bdepth/8);
								buffer.order(ByteOrder.LITTLE_ENDIAN);
								short[] data = (short[])ip.crop().getPixels();
								for(short e : data)buffer.putShort(e);
								System.arraycopy(buffer.array(), 0, bdata, bsize, buffer.array().length);
								bsize += buffer.array().length;
							}else if(bdepth == 32){
								ByteBuffer buffer = ByteBuffer.allocate(b.w_*b.h_*bdepth/8);
								buffer.order(ByteOrder.LITTLE_ENDIAN);
								float[] data = (float[])ip.crop().getPixels();
								for(float e : data)buffer.putFloat(e);
								System.arraycopy(buffer.array(), 0, bdata, bsize, buffer.array().length);
								bsize += buffer.array().length;
							}
						}
						
						String filename = basename + "_Lv" + String.valueOf(l) +
													 "_Ch" + String.valueOf(ch) +
													 "_Fr" + String.valueOf(f) +
													 "_ID" + String.valueOf(i);
						int offset = bytecount;
						int datasize = bdata.length;

						if (filetype == "RAW") {
							int dummy = -1;
							//do nothing
						}
						if (filetype == "JPEG" && bdepth == 8) {
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
								//ImageIO.write(img, format, baos);
								bdata = baos.toByteArray();
								datasize = bdata.length;
							} catch (IOException e) {
								e.printStackTrace();
								return;
							}
						}
						if (filetype == "ZLIB") {
							byte[] tmpdata = new byte[b.w_*b.h_*b.d_*bdepth/8];
							Deflater compresser = new Deflater();
							compresser.setInput(bdata);
							compresser.setLevel(Deflater.DEFAULT_COMPRESSION);
							compresser.setStrategy(Deflater.DEFAULT_STRATEGY);
							compresser.finish();
							datasize = compresser.deflate(tmpdata);
							bdata = tmpdata;
							compresser.end();
						}

						if (bytecount + datasize > sizelimit) {
							BufferedOutputStream fis = null;
							try {
								File file = new File(directory + current_dataname);
								fis = new BufferedOutputStream(new FileOutputStream(file));
								fis.write(packed_data, 0, bytecount);
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
							filecount++;
							current_dataname = base_dataname+"_data"+filecount;
							bytecount = 0;
							offset = 0;
							System.arraycopy(bdata, 0, packed_data, bytecount, datasize);
							bytecount += datasize;
						} else {
							System.arraycopy(bdata, 0, packed_data, bytecount, datasize);
							bytecount += datasize;
						}

						Element filenode = doc.createElement("File");
						filenode.setAttribute("filename", current_dataname);
						filenode.setAttribute("channel", String.valueOf(ch));
						filenode.setAttribute("frame", String.valueOf(f));
						filenode.setAttribute("brickID", String.valueOf(i));
						filenode.setAttribute("offset", String.valueOf(offset));
						filenode.setAttribute("datasize", String.valueOf(datasize));

						fsnode.appendChild(filenode);
	
						curbricknum++;
						IJ.showProgress((double)(curbricknum)/(double)(totalbricknum));

					}
					if (bytecount > 0) {
						BufferedOutputStream fis = null;
						try {
							File file = new File(directory + current_dataname);
							fis = new BufferedOutputStream(new FileOutputStream(file));
							fis.write(packed_data, 0, bytecount);
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
					}
					
				}
			}

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
		
		File newXMLfile = new File(directory + basename + ".vvd");
		writeXML(newXMLfile, doc);
		
		for (int l = 1; l < lv; l++) {
			imp = WindowManager.getImage(lvImgTitle.get(l));
			imp.changes = false;
			imp.close();
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
}
