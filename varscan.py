import urllib
import json
from pyliftover import LiftOver

# http://www.broadinstitute.org/oncotator/mutation/1_226064403_226064403_A_C/
#  in OncoTator JSON


def url_request(chrom, coord, ref, var):
    url = ("http://www.broadinstitute.org/oncotator/mutation/" +
           chrom + "_" +
           coord + "_" +
           coord + "_" +
           ref + "_" +
           var + "/")

    try:
        result = json.loads(
            urllib.request
            .urlopen(url)
            .read()
            .decode('UTF-8')
        )

        if "protein_change" in result:
            return result["protein_change"]
        else:
            raise urllib.error.URLError("Expected JSON not received.")
    except (urllib.error.HTTPError,
            urllib.error.URLError) as e:
        print(e.reason)
        print(url)


def main():
    # do varscan, get ref and var and coords
    ref = "C"
    var = "A"
    hg38_coord = 90088605
    hg38_chrom = "chr15"

    lo = LiftOver("hg38ToHg19.over.chain.gz")
    result = lo.convert_coordinate(hg38_chrom, hg38_coord)

    if result is not None:
        coords_list = result[0]
        protein_change = url_request(coords_list[0],
                                     str(coords_list[1]),
                                     ref,
                                     var)
        print(protein_change)


# returns None or a list
if __name__ == "__main__":
    main()
