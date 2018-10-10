# eoscraper
SCRAPE random numbers generation on EOS

## Short description
A proof-of-concept smart-contract for EOS, allowing users to "book" the random number and pay for it when random number is ready to be reconstructed on client side.


## Description
Existing solutions for random numbers generation in blockchains are based mostly on two main schemes: oracles, providing random numbers and many kinds of "commit-reveal" scheme, when multiple users first commit hashes of generated random numbers in contract, revealing the number itself later. First variant suffers from subverted oracles, the second one - from situation where malicious players "commit", and then decide, to do "reveal" or not, choosing from two variants of resulting random (or more variants, when using sybil attacks). Usage of any blockchain-based parameters (such as blockhash or timestamp) as any part of random seed automatically gives malicious miner ability to choose between N variants of resulting random at the price (N - 1) * block_reward.


### Описание решения
Планируемое решение основано на статье "SCRAPE - Scalable Randomness Attested by Public Entities", основано на протоколе разделения секрета Фиата-Шамира и протокола подбрасывания монеты. Оно обладает следующими свойствами:
 - позволяет раскрыть одно и только одно значение beacon, зафиксированное генерильщиками на этапе фазы commit
 - позволяет раскрыть значение beacon при наличии любого количества честных генерильщиков больше заданного количества, например > 50%
 - позволяет генерировать честный beacon если хотя бы один из генерильщиков честный

Чтобы избавиться от проблемы с нераскрытием участников, был разработан экономический протокол, работающий так:
 - казино "заказывает" рандом, одновременно закладывая в контракт будущую оплату за него
 - казино указывает параметры генерации, регистрирует публичные ключи участников, с которыми хочет работать
 - каждый генерильщик генерирует свой рандом, делит его на shares, шифрует публичными ключами участников и публикует их shares вместе с доказательством, что эти shares были правильно зашифрованы. Вместе с shares генерильщик отправляет депозит, который будет утерян (распределен между честными генерильщиками), если рандом не будет получен за определенный период.
 - генерильщики расшифровывают свои shares. Если кол-во раскрытых shares перевалило за threshold, рандом считается полученным, и его можно вычислить. В этот момент участники протокола получают оплату, притом те участники которые первыми опубликовали расшифрованные шары получают больше прибыли, это мотивирует участников работать быстрее и ускоряет процесс раскрытия.
 - результирующий рандом получается из вскрытых рандомов участников


Дальнейшее ускорение генерации заключается в большом количестве заранее закоммиченных shares.

