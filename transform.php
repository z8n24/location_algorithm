<?php
class transform{
    const PI =3.1415926535898;
    const A = 6378245;
    const EE = 0.006693421622965823;
    // 粗暴判断是否是中国 矩形块。
    public function isInChinaBbox($lon, $lat) :bool{
        return $lon >= 72.004 && $lon <= 137.8347 && $lat >= 0.8293 && $lat <= 55.8271;
    }

    public function transformLat($x , $y ): int {
        $ret = -100 + 2 * $x + 3 * $y + 0.2 * $y * $y + 0.1 * $x * $y + 0.2 * sqrt(abs($x));
        $ret += (20 * sin(6 * $x * self::PI) + 20 * sin(2 * $x * self::PI)) * 2 / 3;
        $ret += (20 * sin($y * self::PI) + 40 * sin($y / 3 * self::PI)) * 2 / 3;
        $ret += (160 * sin($y / 12 * self::PI) + 320 * sin($y * self::PI / 30)) * 2 / 3;
        return $ret;
    }

    public function transformLon($x, $y): int {
        $ret = 300 + $x + 2 * $y + 0.1 * $x * $x + 0.1 * $x * $y + 0.1 * sqrt(abs($x));
          $ret += (20 * sin(6 * $x * self::PI) + 20 * sin(2 * $x * self::PI)) * 2 / 3;
          $ret += (20 * sin($x * self::PI) + 40 * sin($x / 3 * self::PI)) * 2 / 3;
          $ret += (150 * sin($x / 12 * self::PI) + 300 * sin($x / 30 * self::PI)) * 2 / 3;
      return $ret;
    }

    public function delta( $lon, $lat){
        $dLon = $this->transformLon($lon - 105, $lat - 35);
        $dLat = $this->transformLat($lon - 105, $lat - 35);

        $radLat = $lat / 180 * self::PI;
        $magic = sin($radLat);

        $magic = 1 - self::EE * $magic * $magic;
        $sqrtMagic = sqrt($magic);

        $dLon = ($dLon * 180) / (self::A / $sqrtMagic * cos($radLat) * self::PI);
        $dLat = ($dLat * 180) / ((self::A * (1 - self::EE)) / ($magic * $sqrtMagic) * self::PI);

        return [$dLon, $dLat];
    }

    public function WGS84ToGCJ02($lon, $lat){
        if (!$this->isInChinaBbox($lon, $lat)) return [$lon, $lat];

        $d = $this->delta($lon, $lat);
        return [$lon + $d[0], $lat + $d[1]];
    }

    public function GCJ02ToWGS84($lon, $lat){
        if (!$this->isInChinaBbox($lon, $lat)) return [$lon, $lat];

        $wgsLon = $lon;
        $wgsLat =  $lat;

        $tempPoint = $this->WGS84ToGCJ02($wgsLon, $wgsLat);
        $dx = $tempPoint[0] - $lon;
        $dy = $tempPoint[1] - $lat;

        $tempWgsLon = $wgsLon;
        $tempWgsLat = $wgsLat;
        while (abs($dx) > 1e-6 && abs($dy) > 1e-6) {
            $wgsLon -= $dx;
            $wgsLat -= $dy;
            $tempWgsLat = $wgsLat;
            $tempWgsLon = $wgsLon;
            $tempPoint = $this->WGS84ToGCJ02($wgsLon, $wgsLat);
            $dx = $tempPoint[0] - $lon;
            $dy = $tempPoint[1] - $lat;

        }
        return [$wgsLon, $wgsLat];
    }
    //https://gist.github.com/jp1017/71bd0976287ce163c11a7cb963b04dd8
    public function GCJ02ToWGS84_BK($lon, $lat){
        if (!$this->isInChinaBbox($lon, $lat)) return [$lon, $lat];

        $dlat = $this->transformLat($lon - 105.0, $lat - 35.0);
        $dlon = $this->transformLon($lon - 105.0, $lat - 35.0);
        $radlat = $lat / 180.0 * M_PI;
        $magic = sin($radlat);
		$magic = 1 - self::EE * $magic * $magic;
		$sqrtmagic = sqrt($magic);
		$dlat = ($dlat * 180.0) / ((self::A * (1 - self::EE)) / ($magic * $sqrtmagic) * M_PI);
		$dlon = ($dlon * 180.0) / (self::A / $sqrtmagic * cos($radlat) * M_PI);
		$mglat = $lat + $dlat;
		$mglon = $lon+ $dlon;
		return  [$lon * 2 - $mglon, $lat * 2 - $mglat ];
    }

    /**
     * WGS84转GCJ02(火星坐标系)
     *
     * @param lng WGS84坐标系的经度
     * @param lat WGS84坐标系的纬度
     * @return 火星坐标数组
     */
    public function WGS84ToGCJ02_BK($lon, $lat){

        if (!$this->isInChinaBbox($lon, $lat)) return [$lon, $lat];
        $dlat = $this->transformLat($lon - 105.0, $lat - 35.0);
        $dlon = $this->transformLon($lon - 105.0, $lat - 35.0);
        $radlat = $lat / 180.0 * M_PI;
        $magic = sin($radlat);
        $magic = 1 - self::EE * $magic * $magic;
        $sqrtmagic = sqrt($magic);

        $dlat = ($dlat * 180.0) / ((self::A * (1 - self::EE)) / ($magic * $sqrtmagic) * M_PI);
        $dlon = ($dlon * 180.0) / (self::A / $sqrtmagic * cos($radlat) * M_PI);
        $mglat = $lat + $dlat;
        $mglng = $dlon + $lon;
        return [$mglng, $mglat ];
    }
    // 高德的js 装换
    public function gd_distance($lat1 ,$lon1 , $lat2, $lon2){
        $pi80 = M_PI / 180;
        $UQ = 6378137;

        $lat1_ = $lat1 * $pi80;
        $lat2_ = $lat2 * $pi80;
        $uq_  = $UQ * 2;
        $d_ =  $lon2 * $pi80 - $lon1 * $pi80;
        $e = 1- cos($lat2_-$lat1_) + ( 1- cos($d_)) * cos($lat1_) * cos($lat2_) /2;
        $distance = $uq_ * asin(sqrt($e));
        return $distance;
    }

}

