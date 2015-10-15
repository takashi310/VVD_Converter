import ij.*;
import ij.io.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.filter.*;
import ij.plugin.frame.*;

import java.io.*;
import java.nio.*;
import java.util.*;

import javax.xml.parsers.*;
import javax.xml.transform.*;
import javax.xml.transform.dom.*;
import javax.xml.transform.stream.*;

import org.apache.xml.serializer.OutputPropertiesFactory;
import org.w3c.dom.*;

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

public class vvd_export_for_large_data implements PlugIn {
	String basename;
	String directory;
	String tmp_dir;
	String stdir = null;
	String m_prefix = "";
	ArrayList<String> lvImgTitle;
	int bw = 128, bh = 128, bd = 64, lv = 1;
	String filetype = "RAW";
	int jpeg_quality = FileSaver.DEFAULT_JPEG_QUALITY;
	ArrayList<Integer> bwlist = new ArrayList<Integer>();
	ArrayList<Integer> bhlist = new ArrayList<Integer>();
	ArrayList<Integer> bdlist = new ArrayList<Integer>();

	public void run(String arg) {
		if (IJ.versionLessThan("1.49d")) return;
		
		if (!showDialog()) return;
		if (!showDialog2()) return;

		tmp_dir = IJ.getDirectory("Choose a directory to save temp files.");
				
		SaveDialog sd = new SaveDialog("Save as Bricks...", "", "");
		basename = sd.getFileName();
		directory = sd.getDirectory();
		if(basename == null || directory == null) return;
		
		build_bricks();
	}
	
	private boolean showDialog() {
		String[] types = {"RAW", "JPEG"};
		GenericDialog gd = new GenericDialog("Generate Bricks");
		gd.addNumericField("Brick width  (default)",  bw, 0);
		gd.addNumericField("Brick height (default)", bh, 0);
		gd.addNumericField("Brick depth  (default)",  bd, 0);
		gd.addNumericField("Levels", lv, 0);
		gd.addChoice("FileType", types, types[0]);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		bw = (int)gd.getNextNumber();
		bh = (int)gd.getNextNumber();
		bd = (int)gd.getNextNumber();
		lv = (int)gd.getNextNumber();
		filetype = types[gd.getNextChoiceIndex()];
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
			gd.addMessage("Level: " + String.valueOf(i), font);
			gd.addChoice("Image: ", titles, i < tnum ? titles[i] : titles[tnum] );
			gd.addStringField("Brick size (w h d)",  bw + " " + bh + " " + bd);
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

	public void save_splitted_tmp_image( ImagePlus img, ArrayList<Rectangle> rois, String prefix, int z) {

		File tdir = new File(tmp_dir);
		if (tmp_dir == null || !tdir.exists() || rois.isEmpty()) return;

		String basepath = tmp_dir + File.separator;

		basepath += "____tmp_vvd";
		File tdir_f = new File(basepath);
		tdir_f.mkdir();
		stdir = basepath;
		
		if (prefix != null) {
			basepath += File.separator + prefix + String.valueOf(z) + "_";
			m_prefix = prefix;
		}

		int id = 0;
		for (Rectangle r : rois) {
			img.setRoi(r);
			//IJ.log(String.valueOf(r.getX()) + " " + String.valueOf(r.getWidth()));
			ImagePlus tile = new Duplicator().run(img);
			int[] dims = tile.getDimensions();
			String pathname = basepath + String.valueOf(id);
			new FileSaver(tile).saveAsTiff(pathname+".tif");
			//IJ.run(tile, "MHD/MHA compressed ...", "save=" + pathname);
			tile.close();
			id++;
		}

	}
/*
	public ImagePlus load_tmp_image(String name) {
		return IJ.openImage(stdir + File.separator + name + ".tif");
	}
*/

	public void clear_all_tmp_images() {
		File tdir_f = new File(stdir);
		if (!tdir_f.exists()) return;

		File[] files = tdir_f.listFiles();

		for (File f : files) {
			f.delete();
		}
		tdir_f.delete();

		stdir = null;
	}

	public void clear_tmp_images(int z) {
		File tdir_f = new File(stdir);
		if (!tdir_f.exists()) return;

		File[] files = tdir_f.listFiles();

		String target_n = m_prefix + String.valueOf(z) + "_";
		for (File f : files) {
			String fname = f.getName();
			if (fname.startsWith(target_n))
				f.delete();
		}
	}

	public ImagePlus load_tmp_image(int z, int id) {
		return IJ.openImage(stdir + File.separator + m_prefix + String.valueOf(z) + "_" + String.valueOf(id) + ".tif");
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

		long _start = System.currentTimeMillis();

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
			ArrayList<Rectangle> rois = new ArrayList<Rectangle>();
			int mw, mh, md, mw2, mh2, md2;
			double tx0, ty0, tz0, tx1, ty1, tz1;
			double bx0, by0, bz0, bx1, by1, bz1;
			int gridnum = 0;
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

						if (k == 0) rois.add(new Rectangle(x, y, mw2, mh2));
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
			
					long start = System.currentTimeMillis();
					
					Brick b_first = bricks.get(0);
					if(b_first.z_ != 0) IJ.log("warning");
					int st_z = b_first.z_;
					int ed_z = b_first.z_ + b_first.d_;
					for(int s = st_z; s < ed_z; s++) {
						ImagePlus slice = new ImagePlus("slice", stack.getProcessor(imp.getStackIndex(ch+1, s+1, f+1)));
						save_splitted_tmp_image(slice, rois, "z", s);
						slice.close();
					}

					long end = System.currentTimeMillis();
					IJ.log(String.valueOf(end - start)  + " ms");

					for(int i = 0; i < bricks.size(); i++){
						Brick b = bricks.get(i);

						if(ed_z > b.z_ || st_z < b.z_+b.d_){
							if(b.z_ > st_z){
								for(int s = st_z; s < b.z_; s++) clear_tmp_images(s);
								st_z = b.z_;
							}
							else if(b.z_ < st_z){
								IJ.log("warning");
								for(int s = st_z-1; s > b.z_; s--) {
									ImagePlus slice = new ImagePlus("slice", stack.getProcessor(imp.getStackIndex(ch+1, s+1, f+1)));
									save_splitted_tmp_image(slice, rois, "z", s);
									slice.close();
								}
								st_z = b.z_;
							}
							
							if(b.z_+b.d_ > ed_z){
								for(int s = ed_z; s < b.z_+b.d_; s++) {
									ImagePlus slice = new ImagePlus("slice", stack.getProcessor(imp.getStackIndex(ch+1, s+1, f+1)));
									save_splitted_tmp_image(slice, rois, "z", s);
									slice.close();
								}
								ed_z = b.z_+b.d_;
							}
							else if(b.z_+b.d_ < ed_z){
								IJ.log("warning");
								for(int s = ed_z; s > b.z_+b.d_; s--) clear_tmp_images(s);
								ed_z = b.z_+b.d_;
							}
						}
						else {
							IJ.log("warning");
							clear_all_tmp_images();
							st_z = b.z_;
							ed_z = b.z_+b.d_;
							for(int s = st_z; s < ed_z; s++) {
								ImagePlus slice = new ImagePlus("slice", stack.getProcessor(imp.getStackIndex(ch+1, s+1, f+1)));
								save_splitted_tmp_image(slice, rois, "z", s);
								slice.close();
							}
						}

						int gnum = rois.size();
						int bsize = 0;
						byte[] bdata = new byte[b.w_*b.h_*b.d_*(bdepth/8)];
						for (int s = st_z; s < ed_z; s++) {

							ImagePlus bslice = load_tmp_image(s, i % gnum);
							ImageProcessor ip = bslice.getProcessor();
							if(bdepth == 8){
								byte[] data = (byte[])ip.getPixels();
								System.arraycopy(data, 0, bdata, bsize, data.length);
								bsize += data.length;
//								byte[] data = (byte[])ip.crop().getPixels();
//								for(int yy = 0; yy < b.h_; yy++){
//									for(int xx = 0; xx < b.w_; xx++){
//										tsip.putPixel(xx+b.x_, yy+b.y_, data[yy*b.w_+xx]);
//									}
//								}
//								zz++;
							}
							else if(bdepth == 16){
								ByteBuffer buffer = ByteBuffer.allocate(b.w_*b.h_*(bdepth/8));
								buffer.order(ByteOrder.LITTLE_ENDIAN);
								short[] data = (short[])ip.getPixels();
								for(short e : data)buffer.putShort(e);
								System.arraycopy(buffer.array(), 0, bdata, bsize, buffer.array().length);
								bsize += buffer.array().length;
							}else if(bdepth == 32){
								ByteBuffer buffer = ByteBuffer.allocate(b.w_*b.h_*(bdepth/8));
								buffer.order(ByteOrder.LITTLE_ENDIAN);
								float[] data = (float[])ip.getPixels();
								for(float e : data)buffer.putFloat(e);
								System.arraycopy(buffer.array(), 0, bdata, bsize, buffer.array().length);
								bsize += buffer.array().length;
							}
							
							bslice.close();
						}
						
						String filename = basename + "_Lv" + String.valueOf(l) +
													 "_Ch" + String.valueOf(ch) +
													 "_Fr" + String.valueOf(f) +
													 "_ID" + String.valueOf(i);
					
						if (filetype == "RAW") {
							filename += ".raw";
							BufferedOutputStream fis = null;
							try {
								File file = new File(directory + filename);
								fis = new BufferedOutputStream(new FileOutputStream(file));
								fis.write(bdata);
								
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
						if (filetype == "JPEG" && bdepth == 8) {
							filename += ".jpg";
							ByteProcessor bip = new ByteProcessor(b.w_, b.h_*b.d_, bdata);
							ImagePlus biplus = new ImagePlus("", bip);
							FileSaver fs = new FileSaver(biplus);
							fs.setJpegQuality(jpeg_quality);
							fs.saveAsJpeg(directory + filename);
							biplus.close();
						}

						Element filenode = doc.createElement("File");
						filenode.setAttribute("filename", filename);
						filenode.setAttribute("channel", String.valueOf(ch));
						filenode.setAttribute("frame", String.valueOf(f));
						filenode.setAttribute("brickID", String.valueOf(i));
						fsnode.appendChild(filenode);
	
						curbricknum++;
						IJ.showProgress((double)(curbricknum)/(double)(totalbricknum));

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
		
		File newXMLfile = new File(directory + basename + ".xml");
		writeXML(newXMLfile, doc);

		clear_all_tmp_images();

		long _end = System.currentTimeMillis();
		IJ.log("Total Time: " + String.valueOf((_end - _start)/1000.0)  + " sec");

		imp.show();
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
