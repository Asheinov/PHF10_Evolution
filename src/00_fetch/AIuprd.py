import requests
import pandas as pd

def get_aiupred(accession, analysis_type="binding"):

    url = "https://aiupred.elte.hu/rest_api"

    r = requests.get(
        url,
        params={
            "accession": accession,
            "analysis_type": analysis_type
        },
        timeout=60
    )

    r.raise_for_status()

    return r.json()

from pprint import pprint
if __name__ == "__main__":

    data = get_aiupred("Q8WUB8")

    print(data.keys())
    pprint(data["AIUPred"][:5])

    print()

    pprint(data["AIUPred-binding"][:5])