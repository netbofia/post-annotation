

import string
alphabet = list(string.ascii_lowercase)
tebahpla = alphabet[::-1]

decoder = dict(letter for letter in zip(alphabet, tebahpla))

msg = "olleh"

decoded_msg = "".join([decoder[letter] for letter in msg])

print(decoded_msg)