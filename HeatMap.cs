using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Xml.Linq;

namespace HeatMapCs
{
    public class Coordinate
    {
        public double x { get; set; }
        public double y { get; set; }

        public Coordinate(double x, double y)
        {
            this.x = x;
            this.y = y;
        }

        public Coordinate copy()
        {
            return new Coordinate(x, y);
        }

        public override bool Equals(Object obj)
        {
            Coordinate coord = obj as Coordinate;
            if (coord == null)
                return false;

            return x == coord.x && y == coord.y;
        }

        public override int GetHashCode()
        {
            unchecked // Overflow is fine, just wrap
            {
                int hash = 17;
                // Suitable nullity checks etc, of course :)
                hash = hash * 23 + x.GetHashCode();
                hash = hash * 23 + y.GetHashCode();
                return hash;
            }
        }
    }

    public class LatLon : Coordinate
    {
        public double lat
        {
            get { return y; }
            set { y = value; }
        }

        public double lon
        {
            get { return x; }
            set { x = value; }
        }

        public LatLon(double lon, double lat) : base(lon, lat)
        {
            this.lat = lat;
            this.lon = lon;
        }

        public new LatLon copy()
        {
            return new LatLon(lon, lat);
        }
    }

    public class TrackLog
    {
        /// <summary>
        /// for GPX <trkseg> tags
        /// </summary>
        public class Trkseg
        {
            public List<Trkpt> List { get; set; }
        }

        /// <summary>
        /// for GPX <trkpt> tags
        /// </summary>
        public class Trkpt
        {
            public LatLon coords { get; set; }

            public Trkpt(double lat, double lon)
            {
                coords = new LatLon(lon, lat);
            }
        }

        public void _parse(string filename)
        {
            _segments = new List<Trkseg>();

            var doc = XDocument.Load(filename);
            var gpx = XNamespace.Get("http://www.topografix.com/GPX/1/1");
            var tracks = doc.Descendants(gpx + "trk");

            foreach (var trk in tracks)
            {
                var segments = trk.Descendants(gpx + "trkseg");

                foreach (var seg in segments)
                {
                    var points = seg.Descendants(gpx + "trkpt");
                    var segObj = new Trkseg();
                    segObj.List = new List<Trkpt>();

                    foreach (var p in points)
                    {
                        string latS = p.Attribute("lat").Value;
                        string lonS = p.Attribute("lon").Value;

                        var latD = double.Parse(latS, CultureInfo.InvariantCulture);
                        var lonD = double.Parse(lonS, CultureInfo.InvariantCulture);

                        var pObj = new Trkpt(latD, lonD);
                        segObj.List.Add(pObj);
                    }

                    _segments.Add(segObj);
                }
            }
        }

        public string filename { get; set; }
        public List<Trkseg> _segments { get; set; }

        public TrackLog(string filename)
        {
            this.filename = filename;
        }

        public List<Trkseg> segments()
        {
            _parse(filename);
            return _segments;
        }
    }

    public class MercatorProjection
    {
        public const int EARTH_RADIUS = 6378137; //in meters

        private double? _pixels_per_degree;
        private double _pixels_per_radian;
        public double? pixels_per_degree
        {
            get { return _pixels_per_degree; }
            set
            {
                _pixels_per_degree = value;
                _pixels_per_radian = value.Value * (180 / Math.PI);
            }
        }

        public double? meters_per_pixel
        {
            get { return 2 * Math.PI * EARTH_RADIUS / 360 / pixels_per_degree; }
            set { pixels_per_degree = 2 * Math.PI * EARTH_RADIUS / 360 / value; }
        }

        public bool is_scaled
        {
            get { return pixels_per_degree.HasValue; }
        }

        public LatLon project(LatLon coord)
        {
            double x = coord.lon * pixels_per_degree.Value;
            double y = -_pixels_per_radian * Math.Log(Math.Tan((Math.PI / 4 + Math.PI / 360 * coord.lat)));

            return new LatLon(x, y);
        }

        public LatLon inverse_project(Coordinate coord)
        {
            double lat = (360 / Math.PI * Math.Atan(Math.Exp(-coord.y / _pixels_per_radian)) - 90);
            double lon = coord.x / pixels_per_degree.Value;

            return new LatLon(lon, lat);
        }

        public void auto_set_scale(Extent extent_in, double padding, int? width = null, int? height = null)
        {
            //We'll work large to minimize roundoff error.
            double SCALE_FACTOR = 1000000.0;
            pixels_per_degree = SCALE_FACTOR;
            var extent_out = extent_in.map(project);
            padding *= 2; //padding-per-edge -> padding-in-each-dimension

            try
            {
                double? pixels_per_lat = null;

                if (height.HasValue)
                {
                    pixels_per_degree = ((height - padding) / extent_out.size().y * SCALE_FACTOR);
                    pixels_per_lat = pixels_per_degree;
                }

                if (width.HasValue)
                {
                    pixels_per_degree = ((width - padding) / extent_out.size().x * SCALE_FACTOR);
                    if (height.HasValue)
                        pixels_per_degree = Math.Min(pixels_per_degree.Value, pixels_per_lat.Value);
                }
            }
            catch (DivideByZeroException)
            {
                throw new Exception(@"You need at least two data points for auto scaling.
                Try specifying the scale explicitly (or extent height or width).");
            }
        }
    }

    public class Extent
    {
        public LatLon min { get; set; }
        public LatLon max { get; set; }

        public Extent(List<Coordinate> coords = null, List<LineSegment> shapes = null)
        {
            if (coords != null)
            {
                min = new LatLon(coords.Min(c => c.x), coords.Min(c => c.y));
                max = new LatLon(coords.Max(c => c.x), coords.Max(c => c.y));
            }
            else if (shapes != null)
            {
                from_shapes(shapes);
            }
            else
                throw new Exception("Extent must be initialized");
        }

        public void update(Extent other)
        {
            min.x = Math.Min(min.x, other.min.x);
            min.y = Math.Min(min.y, other.min.y);
            max.x = Math.Max(max.x, other.max.x);
            max.y = Math.Max(max.y, other.max.y);
        }

        public void from_bounding_box(Extent other)
        {
            min = other.min.copy();
            max = other.max.copy();
        }

        public void from_shapes(List<LineSegment> shapes)
        {
            foreach (var s in shapes)
                from_bounding_box(s.extent);

            foreach (var s in shapes)
                update(s.extent);
        }

        public Coordinate[] corners()
        {
            return new Coordinate[] { min, max };
        }

        public Coordinate size()
        {
            return new Coordinate(max.x - min.x, max.y - min.y);
        }

        public void grow(int pad)
        {
            min.x -= pad;
            min.y -= pad;
            max.x += pad;
            max.y += pad;
        }

        public void resize(int? width = null, int? height = null)
        {
            if (width.HasValue)
            {
                max.x += (width.Value - size().x) / 2;
                min.x = max.x - width.Value;
            }

            if (height.HasValue)
            {
                max.y += (height.Value - size().y) / 2;
                min.y = max.y - height.Value;
            }
        }

        public bool is_inside(Coordinate coord)
        {
            return (coord.x >= min.x && coord.x <= max.x && coord.y >= min.y && coord.y <= max.y);
        }

        public Extent map(Func<LatLon, Coordinate> func)
        {
            var coords = new List<Coordinate>();
            coords.Add(func(min));
            coords.Add(func(max));

            return new Extent(coords);
        }
    }

    public class Matrix
    {
        public Dictionary<Coordinate, List<double>> _matrix { get; set; }
        public double decay { get; set; }

        public Matrix()
        {
            _matrix = new Dictionary<Coordinate, List<double>>();
        }

        public Matrix(double decay)
        {
            _matrix = new Dictionary<Coordinate, List<double>>();
            this.decay = decay;
        }

        public Extent extent()
        {
            List<Coordinate> coords = _matrix.Keys.ToList();

            return new Extent(coords);
        }

        public void add(Coordinate coord, double val)
        {
            if (!_matrix.ContainsKey(coord) || _matrix[coord] == null)
                _matrix[coord] = new List<double>();

            _matrix[coord].Add(val);
        }

        public Matrix finalized()
        {
            var m = new Matrix();

            foreach (KeyValuePair<Coordinate, List<double>> c in _matrix)
                m._matrix[c.Key] = reduce(decay, c.Value);

            return m;
        }

        public List<double> reduce(double decay, List<double> values)
        {
            double weight = 1.0;
            double total = 0.0;
            values = values.OrderByDescending(d => d).ToList();

            foreach (var value in values)
            {
                total += value * weight;
                weight *= decay;
            }

            var totalList = new List<double>();
            totalList.Add(total);

            return totalList;
        }
    }

    public class LineSegment
    {
        public LatLon start { get; set; }
        public LatLon end { get; set; }
        public double weight { get; set; }
        public double length_squared { get; set; }
        public Extent extent { get; set; }

        public LineSegment(LatLon start, LatLon end, double weight = 1.0)
        {
            this.start = start;
            this.end = end;
            this.weight = weight;
            length_squared = Math.Pow(this.end.x - this.start.x, 2) +
                Math.Pow(this.end.y - this.start.y, 2);

            var coords = new List<Coordinate>();
            coords.Add(start);
            coords.Add(end);
            extent = new Extent(coords);
        }

        public double distance(Coordinate coord)
        {
            double dx = 0;
            double dy = 0;
            double u = 0;

            try
            {
                dx = end.x - start.x;
                dy = end.y - start.y;

                if (length_squared == 0)
                    throw new DivideByZeroException();

                u = ((coord.x - start.x) * dx +
                    (coord.y - start.y) * dy) / length_squared;
                if (u < 0)
                    u = 0;
                else if (u > 1)
                    u = 1;
            }
            catch (DivideByZeroException ex)
            {
                u = 0; // Our line is zero-length.  That's ok.
            }

            dx = start.x + u * dx - coord.x;
            dy = start.y + u * dy - coord.y;

            return Math.Sqrt(dx * dx + dy * dy);
        }

        public void add_heat_to_matrix(Matrix matrix, LinearKernel kernel)
        {
            int fromX = (int)(extent.min.x - kernel.radius);
            int toX = (int)(extent.max.x + kernel.radius + 1);

            int fromY = (int)(extent.min.y - kernel.radius);
            int toY = (int)(extent.max.y + kernel.radius + 1);

            foreach (var x in Main.range(fromX, toX))
            {
                foreach (var y in Main.range(fromY, toY))
                {
                    var coord = new Coordinate(x, y);
                    var dist = distance(coord);
                    double? heat = kernel.heat(dist);

                    if (heat.HasValue && heat != 0)
                        matrix.add(coord, weight * heat.Value);
                }
            }
        }

        public LineSegment map(MercatorProjection projection)
        {
            var x1 = projection.project(start);
            var x2 = projection.project(end);

            return new LineSegment(x1, x2);
        }
    }

    public class LinearKernel
    {
        public int radius { get; set; } //in pixels
        public double radius_float { get; set; }

        public LinearKernel(int radius)
        {
            this.radius = radius;
            radius_float = (double)radius;
        }

        public double heat(double distance)
        {
            if (distance >= radius)
                return 0.0;
            return 1.0 - (distance / radius_float);
        }
    }

    public class ColorMap
    {
        public List<Color> values { get; set; }

        public const string DEFAULT_HSVA_MIN_STR = "000ffff00";
        public const string DEFAULT_HSVA_MAX_STR = "02affffff";

        public static double _str_to_float(string @string, int @base = 16, int maxval = 256)
        {
            int? n = null;
            if (@base == 16)
                n = int.Parse(@string, NumberStyles.HexNumber);
            else
                throw new NotImplementedException("Only hexadecimal numbers");

            return (double)n / maxval;
        }

        public static List<double> str_to_hsva(string @string)
        {
            if (@string.StartsWith("#"))
                @string = @string.Substring(1); //Leading "#" is now optional.

            string s1 = @string.Substring(0, 3);
            string s2 = @string.Substring(3, 3);
            string s3 = @string.Substring(5, 3);
            string s4 = @string.Substring(7, 3);

            var tuple = new List<double>();
            tuple.Add(_str_to_float(s1));
            tuple.Add(_str_to_float(s2));
            tuple.Add(_str_to_float(s3));
            tuple.Add(_str_to_float(s4));

            return tuple;
        }

        public ColorMap(List<double> hsva_min = null, List<double> hsva_max = null, Image image = null, int steps = 256)
        {
            if (image != null)
            {
                bool isRgba = image.PixelFormat == PixelFormat.Format16bppArgb1555 ||
                    image.PixelFormat == PixelFormat.Format32bppArgb ||
                    image.PixelFormat == PixelFormat.Format32bppPArgb ||
                    image.PixelFormat == PixelFormat.Format64bppArgb ||
                    image.PixelFormat == PixelFormat.Format64bppPArgb;

                if (!isRgba)
                    throw new Exception("Gradient image must be RGBA.");

                int num_rows = image.Height;

                using (Bitmap bmp = new Bitmap(image))
                {
                    if (values == null)
                        values = new List<Color>();

                    foreach (var row in Enumerable.Range(0, num_rows))
                    {
                        var c = bmp.GetPixel(0, row);
                        values.Add(c);
                        
                    }

                    values.Reverse();
                }
            }
        }

        public Color get(double floatval)
        {
            return values[(int)(floatval * (values.Count - 1))];
        }
    }

    public class ImageMaker
    {
        public Configuration config { get; set; }
        public Color? background { get; set; }

        public ImageMaker(Configuration config)
        {
            this.config = config;

            if (config.background != null)
                background = config.background;
            else
                background = null;
        }

        public static List<int> _blend_pixels(List<double> a, List<double> b)
        {
            var alpha = a[3] / 255.0;

            var tuple = new List<int>();

            foreach (var i in Enumerable.Range(0, 3))
                tuple[i] = (int)(a[i] * alpha + b[i] * (1 - alpha));

            return tuple;
        }

        public Image make_image(Matrix matrix)
        {
            var extent = config.extent_out;

            if (extent == null)
                extent = matrix.extent();

            int width = config.width.HasValue ? config.width.Value : 1;
            int height = config.height.HasValue ? config.height.Value : 1;

            extent.resize(width - 1, height - 1);

            var size = extent.size();
            size.x = (int)(size.x) + 1;
            size.y = (int)(size.y) + 1;

            var img = new Bitmap((int)size.x, (int)size.y);

            var maxval = matrix._matrix.Max(m => m.Value[0]);

            foreach (KeyValuePair<Coordinate, List<double>> i in matrix._matrix)
            {
                var coord = i.Key;
                double val = i.Value[0];

                if (val == 0.0)
                    continue;

                var x = (int)(coord.x - extent.min.x);
                var y = (int)(coord.y - extent.min.y);

                if (extent.is_inside(coord))
                {
                    var color = config.colormap.get(val / maxval);
                    img.SetPixel(x, y, color);
                }
            }

            return img;
        }
    }

    public class GPXFileReader
    {
        public List<string> filenames { get; set; }

        public GPXFileReader(List<string> filenames)
        {
            this.filenames = filenames;
        }

        public List<LineSegment> iter()
        {
            var segments = new List<LineSegment>();
            foreach (var f in filenames)
            {
                var track = new TrackLog(f);

                foreach (var trkseg in track.segments())
                {
                    if (trkseg.List.Count < 2)
                        continue;

                    foreach (var i in Enumerable.Range(0, trkseg.List.Count - 1))
                    {
                        var p1 = trkseg.List[i];
                        var p2 = trkseg.List[i + 1];

                        segments.Add(new LineSegment(p1.coords, p2.coords));
                    }
                }
            }

            return segments;
        }
    }

    public class Configuration
    {
        public Color? background { get; set; }
        public string output { get; set; }
        public string gradient { get; set; }
        public double decay { get; set; }
        public List<LineSegment> shapes { get; set; }
        public Extent extent_in { get; set; }
        public Extent extent_out { get; set; }
        public ColorMap colormap { get; set; }
        public MercatorProjection projection { get; set; }
        public List<string> files { get; set; }
        public LinearKernel kernel { get; set; }
        public int radius { get; set; }
        public int? width { get; set; }
        public int? height { get; set; }
        public double margin { get; set; }

        public Configuration()
        {
            background = Color.FromArgb(255, 255, 255, 255);
            output = "output.png";
            gradient = "C:\\heatmap\\g.png"; //path to g.png
            decay = 0.95;
            projection = new MercatorProjection();
            colormap = new ColorMap(image: Image.FromFile(gradient));
            radius = 5;
            width = 5000;
            height = 5000;
            kernel = new LinearKernel(radius);
            load_gpx();
        }

        public void load_gpx()
        {
            files = new List<string>();
            string directoryPath = "C:\\heatmap"; //path to input gpx files

            string[] fileEntries = Directory.GetFiles(directoryPath);
            foreach (string fileName in fileEntries)
            {
                if (fileName.EndsWith(".gpx"))
                    files.Add(fileName);
            }
        }

        public void load_files()
        {
            if (files != null)
                shapes = new GPXFileReader(files).iter();
        }

        public void fill_missing()
        {
            var padding = margin + kernel.radius;

            if (extent_in == null)
                extent_in = new Extent(shapes: shapes);

            if (!projection.is_scaled)
                projection.auto_set_scale(extent_in, padding, width, height);
        }
    }

    public static class Log
    {
        public static List<string> log = new List<string>();
    }

    public class Main
    {
        public Image image { get; set; }
        public List<string> log { get; set; }

        public static IEnumerable<int> range(int start, int stop)
        {
            return Enumerable.Range(start, Math.Abs(stop - start));
        }

        public Matrix process_shapes(Configuration config)
        {
            var matrix = new Matrix(config.decay);

            var shapes2 = new List<LineSegment>();

            for (int i = 0; i < config.shapes.Count; i++)
            {
                config.shapes[i] = config.shapes[i].map(config.projection);
                config.shapes[i].add_heat_to_matrix(matrix, config.kernel);
            }

            return matrix;
        }

        public Main()
        {
            var config = new Configuration();
            config.load_files();
            config.fill_missing();
            var matrix = process_shapes(config);
            matrix = matrix.finalized();

            image = new ImageMaker(config).make_image(matrix);
            image.Save(@"C:\heatmap\output.png", ImageFormat.Png); //output path
            log = Log.log;
        }
    }
}
